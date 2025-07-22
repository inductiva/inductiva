
const BANNER_BIG_MESSAGES = {
1: {
    active: true,
    top_text: "Ready to dive in?",
    bot_text: "Get started for free today and earn $5 in credits.",
    button_text: "Get Started",
},
};

const BANNER_SMALL_MESSAGES = {
1: {
    active: true,
    text: "<strong>Start running simulations seamlessly!</strong> You have $5 in <strong>free</strong> credits, no credit card required.",
    button_text: "Sign In",
},
2: {
    active: true,
    text: "<strong>Launch simulations today</strong> with <strong>$5 free credits</strong> - no card needed.",
    button_text: "Get $5 Free",
},
3: {
    active: true,
    text: "<strong>No more queues.</strong> Run simulations instantly with <strong>$5 free.</strong>",
    button_text: "Simulate Now",
},
};

document.addEventListener("DOMContentLoaded", function () {
// Big banner logic
const bigActiveKeys = Object.keys(BANNER_BIG_MESSAGES).filter(
    (key) => BANNER_BIG_MESSAGES[key].active
);
if (bigActiveKeys.length > 0) {
    const selectedBigKey = bigActiveKeys[Math.floor(Math.random() * bigActiveKeys.length)];
    window.banner_text_id = selectedBigKey;
    window.banner_type = "big";
    const selectedBigBanner = BANNER_BIG_MESSAGES[selectedBigKey];

    const headlineElement = document.querySelector("p.headline");
    if (headlineElement) {
    headlineElement.innerHTML = selectedBigBanner.top_text;
    }

    const subtextElement = document.querySelector("p.subtext");
    if (subtextElement) {
    subtextElement.innerHTML = selectedBigBanner.bot_text;
    }

    const bigButtonElement = document.querySelector(".buttons .btn-main");
    if (bigButtonElement) {
    bigButtonElement.innerHTML = selectedBigBanner.button_text;
    }
}

// Small banner logic
const smallActiveKeys = Object.keys(BANNER_SMALL_MESSAGES).filter(
    (key) => BANNER_SMALL_MESSAGES[key].active
);
if (smallActiveKeys.length > 0) {
    const selectedSmallKey = smallActiveKeys[Math.floor(Math.random() * smallActiveKeys.length)];
    window.banner_text_id = selectedSmallKey;
    window.banner_type = "small";
    const selectedSmallBanner = BANNER_SMALL_MESSAGES[selectedSmallKey];

    const ctaTextElement = document.querySelector(".cta-bar .cta-text");
    if (ctaTextElement) {
    ctaTextElement.innerHTML = selectedSmallBanner.text;
    }

    const smallButtonElement = document.querySelector(".cta-bar .cta-button");
    if (smallButtonElement) {
    smallButtonElement.innerHTML = selectedSmallBanner.button_text;
    }
}
});
