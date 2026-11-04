// Vendored and maintained locally for the chart-forge preview loop.
// The original vendor chrome (theme toggle, search box, toolbar buttons,
// footer links and status toasts) was removed: a delivered figure page
// shows the figure, nothing else.

function createHeader(title) {
  return `
        <div class="preview-header">
            <div class="preview-header-info">
                <strong class="preview-header-title">${title}</strong>
            </div>
        </div>
    `;
}
