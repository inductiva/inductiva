window.addEventListener("DOMContentLoaded", function () {
  const themeButton = document.querySelector("nav button[aria-label='Color theme switcher']");
  if (themeButton && !document.getElementById("discord-icon")) {
    const img = document.createElement("img");
    img.src = "_static/Discord-Symbol-Blurple.svg"; // Update if needed
    img.alt = "Discord";
    img.id = "discord-icon";
    img.style.cssText = "height: 24px; margin-right: 8px; cursor: pointer;";
    img.onclick = () => window.open("https://discord.gg/tQ5T3ATt9E", "_blank");

    // Insert before the theme button (i.e., between ⌘K and theme toggle)
    themeButton.parentElement.insertBefore(img, themeButton);
  }
});
