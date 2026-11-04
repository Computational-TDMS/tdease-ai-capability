// Paul Tol colour schemes (SRON technical note; CVD-safe qualitative sets) —
// materialized from pip tol-colors 2.2.0. Use like crameri-*.js:
//   d3.scaleOrdinal().range(TOL_BRIGHT) — or feed d3.interpolateRgbBasis.
// CVD-safe as qualitative sets: bright / vibrant / muted / high_contrast /
// medium_contrast. light is NOT CVD-safe (kept for pale/print contexts).

export const TOL_BRIGHT = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB"];

export const TOL_HIGH_CONTRAST = ["#000000", "#004488", "#BB5566", "#DDAA33", "#FFFFFF"];

export const TOL_MEDIUM_CONTRAST = ["#FFFFFF", "#6699CC", "#004488", "#EECC66", "#997700", "#EE99AA", "#994455", "#000000"];

export const TOL_VIBRANT = ["#EE7733", "#0077BB", "#33BBEE", "#EE3377", "#CC3311", "#009988", "#BBBBBB"];

export const TOL_MUTED = ["#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#44AA99", "#999933", "#AA4499", "#DDDDDD"];

export const TOL_LIGHT = ["#77AADD", "#EE8866", "#EEDD88", "#FFAABB", "#99DDFF", "#44BB99", "#BBCC33", "#AAAA00", "#DDDDDD"];
