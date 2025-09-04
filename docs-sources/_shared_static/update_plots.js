// _static/js/darkmode-plotly.js
document.addEventListener("DOMContentLoaded", () => {
  // Initialize darkMode from localStorage
  let darkMode = localStorage.getItem("theme") || "light";

  // Apply theme to Plotly plots
  function applyPlotlyTheme() {
    if (!window.Plotly) return; // Skip if Plotly not loaded
    const newTemplate = darkMode === "dark" ? "plotly_dark" : "plotly_white";
    document.querySelectorAll(".js-plotly-plot").forEach(plot => {
      Plotly.update(plot, {}, { template: newTemplate });
    });
  }

  // Initial application
  applyPlotlyTheme();

  // Hook up your toggle button
  const btn = document.getElementById("themeToggleButton");
  if (btn) {
    btn.addEventListener("click", () => {
      darkMode = darkMode === "light" ? "dark" : "light";
      localStorage.setItem("theme", darkMode);
      document.documentElement.classList.toggle("dark-mode", darkMode === "dark");
      applyPlotlyTheme();
    });
  }
});
