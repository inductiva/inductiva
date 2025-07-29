
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
        active: false,
        text: "Run your first simulation with <strong>$5 in free credits</strong> — no card needed.",
        button_text: "Claim Free Credits",
    },
    2: {
        active: false,
        text: "<strong>Speed up your simulations</strong> — get <strong>$5 free</strong> to start instantly.",
        button_text: "Start for Free",
    },
    3: {
        active: false,
        text: "<strong>Tired of waiting?</strong> Skip queues and simulate instantly with <strong>$5 free</strong>.",
        button_text: "Simulate Now",
    },
    4: {
        active: true,
        text: "From your laptop to large-scale simulations with five lines of Python.",
        button_text: "Start For Free",
    },
    5: {
        active: true,
        text: "Spin up your personal cluster in minutes and boost your productivity 10x.",
        button_text: "Start For Free",
    },
    6: {
        active: true,
        text: "<strong>Inductiva</strong> runs the simulations, so you can run the science.",
        button_text: "Start For Free",
    },
    7: {
        active: true,
        text: "Leverage the power of high-performance computing in minutes.",
        button_text: "Start For Free",
    },
    8: {
        active: true,
        text: "You can finally simulate as much as you <em>wish</em> you could, not just as much as you can.",
        button_text: "Start For Free",
    },
    9: {
        active: true,
        text: "Turn days of simulation babysitting into minutes of automated, parallel computing.",
        button_text: "Start For Free",
    },
    10: {
        active: true,
        text: "Still juggling simulations between two machines and a prayer? <strong>Inductiva</strong> helps you scale your simulations in minutes.",
        button_text: "Start For Free",
    },
    11: {
        active: true,
        text: "If your laptop is your bottleneck, <strong>Inductiva</strong> is your breakthrough. Instant, on-demand scale with a pay-as-you-go approach.",
        button_text: "Start For Free",
    },
    12: {
        active: true,
        text: "Why wait for cluster time when others are publishing? Run your simulations now and stay ahead.",
        button_text: "Start For Free",
    },
    13: {
        active: true,
        text: "Science moves fast. With <strong>Inductiva</strong>, so can you.",
        button_text: "Start For Free",
    },
    14: {
        active: true,
        text: "Run your simulations only when you need to. No licenses, no subscription traps.",
        button_text: "Start For Free",
    },
    15: {
        active: true,
        text: "You’ve mastered your domain. Let us handle the infrastructure.",
        button_text: "Start For Free",
    },
    16: {
        active: true,
        text: "<strong>Inductiva</strong> offers HPC as it should be: click, run, results.",
        button_text: "Start For Free",
    },
    17: {
        active: true,
        text: "Run massive simulations when your client needs it. Pay nothing when they don’t.",
        button_text: "Start For Free",
    },
    18: {
        active: true,
        text: "Run more simulations. Spend less time configuring.",
        button_text: "Start For Free",
    },
    19: {
        active: true,
        text: "The compute backbone for domain-driven simulation apps.",
        button_text: "Start For Free",
    },
    20: {
        active: true,
        text: "From simulation to dataset model — all in one Python workflow.",
        button_text: "Start For Free",
    },
    21: {
        active: true,
        text: "Forget about queues, job schedulers, or cluster configs. Just import, submit, and simulate.",
        button_text: "Start For Free",
    },
    22: {
        active: true,
        text: "Run simulations like you're at MIT or CERN — from anywhere in the world.",
        button_text: "Start For Free",
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
