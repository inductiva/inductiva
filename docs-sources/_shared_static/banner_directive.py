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
                <div class="text-section">
                    <h1>Ready to dive in?</h1>
                    <p>Get started for free today and earn <span class="highlight">$5 in credits.</span></p>
                    <div class="features">
                    <span>• No credit card</span>
                    <span>• Free forever</span>
                    </div>
                </div>
                <div class="button-section">
                    <button class="get-started-btn">Get Started</button>
                </div>
                </div>
            </div>
            '''
        return [nodes.raw('', html, format='html')]
