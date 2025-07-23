window.addEventListener('scroll', () => {
  const buffer = 100; // px from the bottom
  const scrollY = window.scrollY || window.pageYOffset;
  const viewportHeight = window.innerHeight;
  const documentHeight = document.documentElement.scrollHeight;

  const nearBottom = (scrollY + viewportHeight) >= (documentHeight - buffer);
  document.body.classList.toggle('at-bottom', nearBottom);
});