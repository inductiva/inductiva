const googleanalytics_id = window.env?.GTAG_WEBSITE || "";


(function(w,d,s,l,i){w[l]=w[l]||[];w[l].push({'gtm.start':
new Date().getTime(),event:'gtm.js'});var f=d.getElementsByTagName(s)[0],
j=d.createElement(s),dl=l!='dataLayer'?'&l='+l:'';j.async=true;j.src=
'https://www.googletagmanager.com/gtm.js?id='+i+dl;f.parentNode.insertBefore(j,f);
})(window,document,'script','dataLayer',googleanalytics_id);

//GTM variables
//Get main frame url
const params = new URLSearchParams(window.top.location.search);
window.url_referrer = params.get('guides_cta_origin') || 'None';