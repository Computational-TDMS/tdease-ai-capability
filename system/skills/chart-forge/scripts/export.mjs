#!/usr/bin/env node
/* Validate self-contained chart-forge HTML in Chromium, then export PNG and,
 * for SVG figures, SVG. Any page, readiness, network, or rendering error is a
 * failed export.
 */
import { createRequire } from "node:module";
import { mkdirSync, writeFileSync, existsSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

const FORGE = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const localRequire = createRequire(import.meta.url);
const cwdRequire = createRequire(resolve(process.cwd(), "package.json"));
const HELP = `Usage:
  node export.mjs --in figure.html [--in other.html] --out export
                  [--scale 2] [--width 1054] [--height 760] [--timeout 20000]`;

function loadPlaywright() {
  const explicit = process.env.CHART_FORGE_PLAYWRIGHT;
  const vendored = resolve(FORGE, "assets/node_modules/playwright");
  for (const path of [explicit, vendored]) {
    if (path && existsSync(path)) return localRequire(path);
  }
  for (const resolver of [localRequire, cwdRequire]) {
    for (const name of ["playwright", "@playwright/test"]) {
      try { return resolver(name); } catch { /* try the next supported location */ }
    }
  }
  throw new Error(
    "Playwright not found. Run from a workspace that provides playwright or " +
    "set CHART_FORGE_PLAYWRIGHT to its package directory.");
}

function argParse(argv) {
  const args = { in: [], out: "export", scale: 2, width: 1054, height: 760, timeout: 20000 };
  const value = (index, flag) => {
    if (index + 1 >= argv.length) throw new Error(`${flag} requires a value`);
    return argv[index + 1];
  };
  for (let i = 0; i < argv.length; i++) {
    const key = argv[i];
    if (key === "--in") args.in.push(resolve(value(i++, key)));
    else if (key === "--out") args.out = resolve(value(i++, key));
    else if (key === "--scale") args.scale = Number(value(i++, key));
    else if (key === "--width") args.width = Number(value(i++, key));
    else if (key === "--height") args.height = Number(value(i++, key));
    else if (key === "--timeout") args.timeout = Number(value(i++, key));
    else if (key === "--help" || key === "-h") { console.log(HELP); process.exit(0); }
    else throw new Error(`unknown argument: ${key}`);
  }
  if (!args.in.length) throw new Error("at least one --in is required");
  for (const key of ["scale", "width", "height", "timeout"]) {
    if (!Number.isFinite(args[key]) || args[key] <= 0) throw new Error(`--${key} must be positive`);
  }
  for (const input of args.in) if (!existsSync(input)) throw new Error(`input not found: ${input}`);
  return args;
}

const args = argParse(process.argv.slice(2));
const { chromium } = loadPlaywright();
mkdirSync(args.out, { recursive: true });

const browser = await chromium.launch();
const results = [];

try {
  for (const input of args.in) {
    const stem = input.replace(/\.html$/i, "").split(/[\\/]/).pop();
    let page;
    try {
      page = await browser.newPage({
        viewport: { width: args.width, height: args.height },
        deviceScaleFactor: args.scale,
      });
      const pageErrors = [];
      const networkRequests = [];
      page.on("pageerror", error => pageErrors.push(String(error).split("\n")[0]));
      page.on("request", request => {
        if (/^https?:/i.test(request.url())) networkRequests.push(request.url());
      });

    await page.goto(pathToFileURL(input).href, { waitUntil: "load", timeout: args.timeout });
    await page.waitForFunction(
      "window.__chartReady === true || Boolean(window.__chartError)",
      null,
      { timeout: args.timeout },
    );
    if (pageErrors.length) throw new Error(`page error: ${pageErrors[0]}`);
    if (networkRequests.length) throw new Error(`external request: ${networkRequests[0]}`);

    const visual = await page.evaluate(() => {
      const graphic = document.querySelector("#chart svg, #chart canvas");
      if (!graphic) return null;
      const box = graphic.getBoundingClientRect();
      return {
        tag: graphic.tagName.toLowerCase(),
        width: box.width,
        height: box.height,
        marks: graphic.tagName.toLowerCase() === "svg" ? graphic.querySelectorAll("*").length : null,
        error: window.__chartError || null,
      };
    });
    if (!visual) throw new Error("no SVG or canvas in #chart");
    if (visual.error) throw new Error(`chart error: ${visual.error}`);
    if (visual.width < 20 || visual.height < 20) {
      throw new Error(`graphic has invalid bounds ${visual.width}x${visual.height}`);
    }
    if (visual.tag === "svg" && visual.marks < 2) throw new Error("SVG has no rendered marks");

    await page.addStyleTag({ content: "#export-bar{display:none!important}" });
    const png = resolve(args.out, `${stem}.png`);
    const app = await page.$("#app");
    if (!app) throw new Error("missing #app root");

    const inner = await page.$("#chart svg, #chart canvas");
    const innerBox = inner && await inner.boundingBox();
    if (innerBox) {
      await page.addStyleTag({ content:
        `#app{width:fit-content}header h1,header p,#provenance{width:${Math.ceil(innerBox.width)}px}` });
      await page.waitForTimeout(80);
    }
    const appBox = await app.boundingBox();
    if (!appBox) throw new Error("#app has no layout box");
    await page.setViewportSize({
      width: Math.max(1, Math.ceil(appBox.width)),
      height: Math.max(1, Math.ceil(appBox.height)),
    });
    await page.waitForTimeout(100);
    await app.screenshot({ path: png });
    results.push(`PNG ${png}`);

    const hasSVG = await page.evaluate(
      () => Boolean(document.querySelector("#chart svg")) &&
        !document.querySelector("#chart canvas") &&
        typeof window.chartForge?.exportSVG === "function");
    if (hasSVG) {
      const svg = await page.evaluate(() => window.chartForge.exportSVG());
      const svgPath = resolve(args.out, `${stem}.svg`);
      writeFileSync(svgPath, svg, "utf8");
      results.push(`SVG ${svgPath}`);
    } else {
      results.push(`SVG skipped ${stem} (canvas/hybrid/webgl render)`);
    }
    } catch (error) {
      results.push(`FAIL ${stem}: ${String(error).split("\n")[0]}`);
    } finally {
      if (page) await page.close();
    }
  }
} finally {
  await browser.close();
}

console.log(results.join("\n"));
if (results.some(result => result.startsWith("FAIL"))) process.exit(1);
