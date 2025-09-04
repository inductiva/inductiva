// _static/js/darkmode-plotly.js
document.addEventListener("DOMContentLoaded", () => {
  // Initialize darkMode from localStorage
  let darkMode = localStorage.getItem("theme") || "light";

  // Apply theme to Plotly plots
  function applyPlotlyTheme() {
    if (!window.Plotly) return; // Skip if Plotly not loaded
    console.log("Im here")
    const newTemplate = darkMode === "dark" ? "plotly_dark" : "plotly_white";
    document.querySelectorAll(".js-plotly-plot").forEach(plot => {
      Plotly.update(plot, {}, { template: newTemplate });
    });
  }

  // Initial application
  applyPlotlyTheme();

  // Hook up your toggle button
  const btn = document.querySelector('button[aria-label="Color theme switcher"]');
  if (btn) {
    btn.addEventListener("click", () => {
        console.log("Clicked theme button");
      darkMode = darkMode === "light" ? "dark" : "light";
      console.log("New theme:", darkMode);
      localStorage.setItem("theme", darkMode);
      document.documentElement.classList.toggle("dark-mode", darkMode === "dark");
      applyPlotlyTheme();
    });
  }
});
