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

    BANNER_MESSAGES = {
        1:{
            "active": True,
            "top_text":"Ready to dive in?",
            "bot_text":"Get started for free today and earn $5 in credits.",
            "button_text":"Get Started"
        },
        2:{
            "active": True,
            "top_text":"Join us!",
            "bot_text":"Register for free today and earn $5 in credits.",
            "button_text":"Get Started For Free!"
        }
    }

    def run(self):
        origin = self.options.get('origin', '')

        #Get only active messages
        active_keys = [key for key, val in self.BANNER_MESSAGES.items() if val.get("active")]

        selected_key = random.choice(active_keys)

        # pylint: disable=line-too-long
        html = f'''
            <div class="banner">
                <div class="banner-content">
                <div class="text">
                    <p class="headline">{self.BANNER_MESSAGES[selected_key]["top_text"]}</p>
                    <p class="subtext">{self.BANNER_MESSAGES[selected_key]["bot_text"]}</p>
                </div>
                <div class="buttons">
                    <button onclick="window.open('https://console.inductiva.ai/?utm_source=guide_{origin}&utm_medium=button&utm_campaign=signup', '_blank')" target="_blank" class="btn primary" id="login-btn-big" >
                    <span class="btn-main">{self.BANNER_MESSAGES[selected_key]["button_text"]}</span>
                    </button>
                </div>
                </div>
            </div>
            <script>
                window.banner_text_id={selected_key};
                window.banner_type="big";
            </script>
            '''
        return [nodes.raw('', html, format='html')]
