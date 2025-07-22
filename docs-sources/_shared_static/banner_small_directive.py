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

    BANNER_MESSAGES = {
        1:{
            "active": True,
            "text":"<strong>Start running simulations seamlessly!</strong> You have $5 in <strong>free</strong> credits, no credit card required.",
            "button_text":"Sign In"
        },
        2:{
            "active": True,
            "text":"<strong>Launch simulations today</strong> with <strong>$5 free credits</strong> - no card needed.",
            "button_text":"Get $5 Free"
        },
        3:{
            "active": True,
            "text":"<strong>No more queues.</strong> Run simulations instantly with <strong>$5 free.</strong>",
            "button_text":"Simulate Now"
        },
    }

    def run(self):
        origin = self.options.get('origin', '')

        #Get only active messages
        active_keys = [key for key, val in self.BANNER_MESSAGES.items() if val.get("active")]

        selected_key = random.choice(active_keys)

        # pylint: disable=line-too-long
        html = f'''
            <div class="cta-bar">
                <div class="cta-text">
                    {self.BANNER_MESSAGES[selected_key]["text"]}
                </div>
                <button  onclick="window.open('https://console.inductiva.ai/?utm_source=guide_{origin}&utm_medium=button&utm_campaign=signup', '_blank')" target="_blank" class="cta-button" id="login-btn-small">{self.BANNER_MESSAGES[selected_key]["button_text"]}</button>
            </div>
            <script>
                window.banner_text_id={selected_key};
                window.banner_type="big";
            </script>
            '''
        return [nodes.raw('', html, format='html')]
