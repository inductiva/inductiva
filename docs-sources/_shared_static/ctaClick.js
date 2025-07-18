

function handleCtaClick(origin) {

    window.dataLayer.push({
        event: 'login_button_click_inside_documentation',
        origin: origin,
        eventCallback: function () {
          callbackCalled = true;
          console.log(`Event callback triggered for origin '${origin}'`);
        }
      });
  
}