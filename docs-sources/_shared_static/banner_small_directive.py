"""Python directive for sphinx banner_small"""
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
                    <strong>Start running simulations seamlessly!</strong> You have $5 in <strong>free</strong> credits, no credit card required.
                </div>
                <button  onclick="window.open('https://console.inductiva.ai/?utm_source=guide_{origin}&utm_medium=button&utm_campaign=signup', '_blank')" target="_blank" class="cta-button">Sign In</button>
            </div>
            '''
        return [nodes.raw('', html, format='html')]
