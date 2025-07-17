window.addEventListener("DOMContentLoaded", function () {
  const themeButton = document.querySelector("nav button[aria-label='Color theme switcher']");

  const pathParts = window.location.pathname.split('/');
  console.log(window.location.pathname)
  console.log(pathParts)
  const simulatorIndex = pathParts.indexOf('guides') + 1;
  const simulator = pathParts[simulatorIndex];

  // Construct the static path
  const staticPath = `/builds/${simulator}/_static/`;

  if (themeButton && !document.getElementById("discord-icon")) {
    const img = document.createElement("img");
    img.src = staticPath+"Discord-Symbol-Blurple.svg";
    img.alt = "Discord";
    img.id = "discord-icon";
    img.style.cssText = `
      height: 24px;
      margin-right: 8px;
      cursor: pointer;
      transition: filter 0.2s ease;
    `;
    img.onclick = () => window.open("https://discord.gg/tQ5T3ATt9E", "_blank");

    img.addEventListener("mouseover", () => {
      img.style.filter = "brightness(1.3)";
    });

    img.addEventListener("mouseout", () => {
      img.style.filter = "brightness(1)";
    });

    themeButton.parentElement.insertBefore(img, themeButton);
  }
});
