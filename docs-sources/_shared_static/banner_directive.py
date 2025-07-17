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
            f'''<div class="banner">
                <h2>FROM: {origin}. Start running simulations seamlessly!</h2>
                <p>You have <strong>$5</strong> in free credits, no credit card required.</p>
                <button class="sign-in-btn">Sign In</button>
            </div>'''
        )
        return [nodes.raw('', html, format='html')]
