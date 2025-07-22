"""Python directive for sphinx banner_small"""
import random
from docutils import nodes
from docutils.parsers.rst import Directive


class BannerSmallDirective(Directive):
    """Python directive for sphinx banner_small"""
    required_arguments = 0
    option_spec = {
        'origin': str,
    }

    def run(self):
        origin = self.options.get('origin', '')

        # pylint: disable=line-too-long
        html = f'''
            <div class="cta-bar">
                <div class="cta-text">
                    ...
                </div>
                <button  onclick="window.open('https://console.inductiva.ai/?utm_source=guide_{origin}&utm_medium=button&utm_campaign=signup', '_blank')" target="_blank" class="cta-button" id="login-btn-small">...</button>
            </div>
            <script>
                const activeKeys = Object.keys(window.BANNER_SMALL_MESSAGES).filter(
                    (key) => BANNER_SMALL_MESSAGES[key].active
                );
                window.banner_text_id = selectedKey = activeKeys[Math.floor(Math.random() * activeKeys.length)];
                window.banner_type="small";
                const selectedBanner = window.BANNER_MESSAGES[window.banner_text_id];
                
                // Replace the text
                const ctaTextElement = document.querySelector(".cta-bar .cta-text");
                if (ctaTextElement) {{
                    ctaTextElement.innerHTML = selectedBanner.text;
                }}
                const buttonElement = document.querySelector(".cta-bar .cta-button");
                if (buttonElement) {{
                    buttonElement.innerHTML = selectedBanner.button_text;
                }}
            </script>
            '''
        return [nodes.raw('', html, format='html')]
