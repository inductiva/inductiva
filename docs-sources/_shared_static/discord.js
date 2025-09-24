window.addEventListener("DOMContentLoaded", function () {
  const themeButton = document.querySelector("nav button[aria-label='Color theme switcher']");
  // get image path based on current simulator to avoid
  // depth problems. use /builds/simulator/_static instead of
  // _/static or ../_static or ../../_static
  const pathParts = window.location.pathname.split('/');
  const simulator = pathParts[2];

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
    img.onclick = () => window.open("https://discord.gg/rFkHxVmAbu", "_blank");

    img.addEventListener("mouseover", () => {
      img.style.filter = "brightness(1.3)";
    });

    img.addEventListener("mouseout", () => {
      img.style.filter = "brightness(1)";
    });

    themeButton.parentElement.insertBefore(img, themeButton);
  }
});
