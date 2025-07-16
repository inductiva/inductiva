window.addEventListener("DOMContentLoaded", function () {
  const searchBar = document.querySelector("input[type='search'], .search-input");
  if (searchBar && !document.getElementById("custom-search-icon")) {
    const img = document.createElement("img");
    img.src = "_static/Discord-Symbol-Blurple.svg";  // Update with your actual image path
    img.alt = "Discord";
    img.id = "discord-icon";
    img.style.cssText = "height: 24px; margin-left: 8px; cursor: pointer;";
    img.onclick = () => window.open("https://discord.gg/tQ5T3ATt9E", "_blank");

    searchBar.parentElement.appendChild(img);
  }
});
