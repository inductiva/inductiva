document.addEventListener("DOMContentLoaded", () => {
  // Initial theme from localStorage or default
  let darkMode = localStorage.getItem("theme") || "dark";

  // Function to update SVG backgrounds
  function updateSVGBackgrounds() {
    const svgs = document.querySelectorAll("svg.main-svg");
    const bgColor = darkMode === "dark" ? "#ffffff7a" : "#ffffff";
    svgs.forEach(svg => {
      svg.style.background = bgColor;
    });
  }

  // Apply initial background
  updateSVGBackgrounds();

  // Hook theme toggle button (assumes aria-label="Color theme switcher")
  const btn = document.querySelector('button[aria-label="Color theme switcher"]');
  if (btn) {
    btn.addEventListener("click", () => {
      // Toggle theme
      darkMode = darkMode === "light" ? "dark" : "light";
      localStorage.setItem("theme", darkMode);

      // Optional: toggle class on <html> for site styles
      document.documentElement.classList.toggle("dark-mode", darkMode === "dark");

      // Update all SVG backgrounds
      updateSVGBackgrounds();
    });
  }
});
