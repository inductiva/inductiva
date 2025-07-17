from docutils import nodes
from docutils.parsers.rst import Directive

class BannerDirective(Directive):
    required_arguments = 0
    option_spec = {
        'origin': str,
    }

    def run(self):
        origin = self.options.get("origin", "")
        html = (
            f'''
            <div class="banner">
                <div class="banner-content">
                <div class="text">
                    <p class="headline">Ready to dive in?</p>
                    <p class="subtext">Start your free trial today and earn $5 in free credits.</p>
                </div>
                <div class="buttons">
                    <button onclick="window.open('https://console.inductiva.ai/?utm_source=guide_{origin}&utm_medium=button&utm_campaign=signup', '_blank')" target="_blank" class="btn primary">
                    <span class="btn-main">Get Started</span>
                    <span class="btn-sub">Free</span>
                    </button>
                </div>
                </div>
            </div>
            '''
        )
        return [nodes.raw('', html, format='html')]
