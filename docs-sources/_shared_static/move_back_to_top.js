window.addEventListener('scroll', () => {
  //Only does it for screens smaller then 1280
  //If the screen is smaller than 1280 we lose the right side
  const isMobile = window.innerWidth <= 1279;
  if (!isMobile) {
    document.body.classList.remove('at-bottom');
    return;
  }

  const buffer = 200;
  const scrollY = window.scrollY || window.pageYOffset;
  const viewportHeight = window.innerHeight;
  const documentHeight = document.documentElement.scrollHeight;

  const nearBottom = (scrollY + viewportHeight) >= (documentHeight - buffer);
  document.body.classList.toggle('at-bottom', nearBottom);
});