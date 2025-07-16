window.addEventListener("DOMContentLoaded", function () {
  const themeButton = document.querySelector("nav button[aria-label='Color theme switcher']");
  if (themeButton && !document.getElementById("discord-icon")) {
    // Create Discord icon
    const img = document.createElement("img");
    img.src = "_static/Discord-Symbol-Blurple.svg";
    img.alt = "Discord";
    img.id = "discord-icon";
    img.style.cssText = `
      height: 24px;
      margin-right: 8px;
      cursor: pointer;
      transition: filter 0.2s ease;
      position: relative;
    `;
    img.onclick = () => window.open("https://discord.gg/tQ5T3ATt9E", "_blank");

    // Create tooltip
    const tooltip = document.createElement("div");
    tooltip.textContent = "Join our Discord to get faster support and chat with us!";
    tooltip.style.cssText = `
      position: absolute;
      bottom: 130%;
      left: 50%;
      transform: translateX(-50%);
      background: #333;
      color: #fff;
      padding: 6px 10px;
      border-radius: 4px;
      font-size: 12px;
      white-space: nowrap;
      opacity: 0;
      pointer-events: none;
      transition: opacity 0.3s ease;
      z-index: 1000;
    `;

    // Create a wrapper to position tooltip relative to icon
    const wrapper = document.createElement("div");
    wrapper.style.cssText = "position: relative; display: inline-block;";
    wrapper.appendChild(img);
    wrapper.appendChild(tooltip);

    // Add hover behavior with delay
    let showTimeout;
    img.addEventListener("mouseenter", () => {
      showTimeout = setTimeout(() => {
        tooltip.style.opacity = "1";
      }, 100);
    });
    img.addEventListener("mouseleave", () => {
      clearTimeout(showTimeout);
      tooltip.style.opacity = "0";
    });

    // Insert icon + tooltip before the theme toggle button
    themeButton.parentElement.insertBefore(wrapper, themeButton);
  }
});
