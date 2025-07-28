
const BANNER_BIG_MESSAGES = {
    1: {
        active: true,
        top_text: "Ready to dive in?",
        bot_text: "Get started for free today and earn $5 in credits.",
        button_text: "Get Started ✨",
    },
};

const BANNER_SMALL_MESSAGES = {
    1: {
        active: true,
        text: "Run your first simulation with <strong>$5 in free credits</strong> — no card needed.",
        button_text: "Claim Free Credits",
    },
    2: {
        active: true,
        text: "<strong>Speed up your simulations</strong> — get <strong>$5 free</strong> to start instantly.",
        button_text: "Start for Free",
    },
    3: {
        active: true,
        text: "<strong>Tired of waiting?</strong> Skip queues and simulate instantly with <strong>$5 free</strong>.",
        button_text: "Simulate Now",
    },
};

document.addEventListener("DOMContentLoaded", function () {
    // Add simulator name to global var
    const sim_name = window.location.pathname.split('/');
    window.simulator_name = sim_name[2];

    // Big banner logic
    const bigActiveKeys = Object.keys(BANNER_BIG_MESSAGES).filter(
        (key) => BANNER_BIG_MESSAGES[key].active
    );
    if (bigActiveKeys.length > 0) {
        const selectedBigKey = bigActiveKeys[Math.floor(Math.random() * bigActiveKeys.length)];
        window.banner_text_id = selectedBigKey;
        window.banner_type = "big";
        const selectedBigBanner = BANNER_BIG_MESSAGES[selectedBigKey];

        document.querySelectorAll("p.headline").forEach(el => {
            el.innerHTML = selectedBigBanner.top_text;
        });

        document.querySelectorAll("p.subtext").forEach(el => {
            el.innerHTML = selectedBigBanner.bot_text;
        });

        document.querySelectorAll(".buttons .btn-main").forEach(el => {
            el.innerHTML = selectedBigBanner.button_text;
        });
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

        document.querySelectorAll(".cta-bar .cta-text").forEach(el => {
            el.innerHTML = selectedSmallBanner.text;
        });

        document.querySelectorAll(".cta-bar .cta-button").forEach(el => {
            el.innerHTML = selectedSmallBanner.button_text;
        });
    }
});
