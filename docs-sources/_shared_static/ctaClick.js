function trackCtaClick(origin, callback) {
  // Ensure gtag is initialized only once
  if (!window.gtagInitialized) {
    console.log('[trackCtaClick] Initializing gtag...');
    
    // Load GA4 script dynamically
    const gtagScript = document.createElement('script');
    gtagScript.async = true;
    gtagScript.src = "https://www.googletagmanager.com/gtag/js?id=G-MJM76GDYEX&cx=c&gtm=45He57g1v9204454354za200&tag_exp=101509157~103116026~103200004~103233427~103351869~103351871~104684208~104684211~104732253~104732255~104908321~104908323~104921715~104921717~104964065~104964067~104967141~104967143";

    gtagScript.onload = () => {
      console.log('[trackCtaClick] GA4 script loaded.');

      window.dataLayer = window.dataLayer || [];
      window.gtag = function () { dataLayer.push(arguments); };

      gtag('js', new Date());

      // ✅ CORRECT GA4 Measurement ID (not GTM)
      gtag('config', 'G-MJM76GDYEX', { debug_mode: true });

      window.gtagInitialized = true;

      // After initialization, send the event
      sendCtaEvent(origin, callback);
    };

    gtagScript.onerror = () => {
      console.error('[trackCtaClick] Failed to load GA4 script.');
      if (callback) callback(false);
    };

    document.head.appendChild(gtagScript);
  } else {
    // Already initialized, just send the event
    sendCtaEvent(origin, callback);
  }
}

// Helper function to send the CTA click event
function sendCtaEvent(origin, callback) {
  if (typeof gtag !== 'function') {
    console.warn('[trackCtaClick] gtag is not defined. Event not sent.');
    if (callback) callback(false);
    return;
  }

  console.log(`[trackCtaClick] Sending event 'cta_click' with origin: ${origin}`);

  gtag('event', 'cta_click', {
    event_category: 'CTA',
    event_label: origin,
    origin: origin
  }, {
    event_callback: function () {
      console.log(`[trackCtaClick] Event 'cta_click' with origin '${origin}' sent successfully.`);
      if (callback) callback(true);
    }
  });
}



function handleCtaClick(origin) {
  trackCtaClick(origin, function(success) {
    const baseUrl = 'https://console.inductiva.ai/?utm_source=guide_' + encodeURIComponent(origin) +
                    '&utm_medium=button&utm_campaign=signup';

    if (success) {
      console.log('[handleCtaClick] Opening signup URL after tracking.');
    } else {
      console.warn('[handleCtaClick] Tracking failed or not initialized — opening anyway.');
    }

    // Open the URL after tracking (or immediately on failure)
    window.open(baseUrl, '_blank');
  });
}