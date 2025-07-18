function trackAndRedirect(origin) {
    // Capture the page as the origin (e.g. /getting-started)
    origin = origin || 'unknown';

    // Compose final redirect URL with UTM tracking
    const url = `https://console.inductiva.ai/`;

    // Track the event in Google Analytics (GA4)
    gtag('event', 'cta_click', {
      origin_page: origin,
      event_category: 'engagement',
      event_label: 'Get Started Button',
      event_callback: function() {
        // Open in new tab after GA sends the event
        window.open(url, '_blank');
      },
      event_timeout: 800
    });
  }