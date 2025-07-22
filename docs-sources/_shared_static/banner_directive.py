"""Python directive for sphinx banner_small"""
import random
from docutils import nodes
from docutils.parsers.rst import Directive


class BannerDirective(Directive):
    """Python directive for sphinx banner_small"""
    required_arguments = 0
    option_spec = {
        'origin': str,
    }

    def run(self):
        origin = self.options.get('origin', '')

        # pylint: disable=line-too-long
        html = f'''
            <div class="banner">
                <div class="banner-content">
                <div class="text">
                    <p class="headline">...</p>
                    <p class="subtext">...</p>
                </div>
                <div class="buttons">
                    <button onclick="window.open('https://console.inductiva.ai/?utm_source=guide_{origin}&utm_medium=button&utm_campaign=signup', '_blank')" target="_blank" class="btn primary" id="login-btn-big" >
                    <span class="btn-main">...</span>
                    </button>
                </div>
                </div>
            </div>
            <script>
                const activeKeys = Object.keys(window.BANNER_BIG_MESSAGES).filter(
                    (key) => BANNER_BIG_MESSAGES[key].active
                );
                window.banner_text_id = selectedKey = activeKeys[Math.floor(Math.random() * activeKeys.length)];
                window.banner_type="big";
                const selectedBanner = window.BANNER_MESSAGES[window.banner_text_id];
                
                // Replace the text
                const headlineElement = document.querySelector("p.headline");
                if (headlineElement) {{
                    headlineElement.innerHTML = selectedBanner.top_text;
                }}
                const subtextElement = document.querySelector("p.subtext");
                if (subtextElement) {{
                    subtextElement.innerHTML = selectedBanner.bot_text;
                }}
                const buttonElement = document.querySelector(".buttons .btn-main");
                if (buttonElement) {{
                    buttonElement.innerHTML = selectedBanner.button_text;
                }}
            </script>
            '''
        return [nodes.raw('', html, format='html')]
