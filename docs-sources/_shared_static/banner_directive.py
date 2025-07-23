"""Python directive for sphinx banner_small"""
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
                    <span class="btn-main">...</span> <span class="sparkle">✨</span>
                    </button>
                </div>
                </div>
            </div>
            '''
        return [nodes.raw('', html, format='html')]
