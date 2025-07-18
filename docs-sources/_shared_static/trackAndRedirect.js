function trackAndRedirect(origin) {
  
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());
  gtag('config', 'G-MJM76GDYEX'); // Your GA4 Measurement ID

  gtag('set', 'debug_mode', true);

  // Capture the page as the origin (e.g. /getting-started)
  origin = origin || 'unknown';
  console.log(`[DEBUG] trackAndRedirect called with origin: "${origin}"`);

  // Compose final redirect URL with UTM tracking
  const url = `https://console.inductiva.ai/`;
  console.log(`[DEBUG] Redirect URL set to: ${url}`);

  // Track the event in Google Analytics (GA4)
  gtag('event', 'cta_click', {
    origin_page: origin,
    event_category: 'engagement',
    event_label: 'Get Started Button',
    event_callback: function() {
      console.log('[DEBUG] GA event_callback fired — opening new tab now.');
      // Open in new tab after GA sends the event
      window.open(url, '_blank');
    },
    event_timeout: 800
  });

  console.log('[DEBUG] gtag event sent, waiting for callback or timeout...');
}

gtag('event', 'test_event', {'app_name': 'myAppName','screen_name': 'Home'});