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
                    <p class="headline">Ready to dive in?</p>
                    <p class="subtext">Get started for free today and earn $5 in credits.</p>
                </div>
                <div class="buttons">
                    <button onclick="handleCtaClick('{origin}')" target="_blank" class="btn primary cta-button">
                    <span class="btn-main">Get Started</span>
                    </button>
                </div>
                </div>
            </div>
            '''
        return [nodes.raw('', html, format='html')]
