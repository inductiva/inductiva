"""Python directive for sphinx banner_small"""
from docutils import nodes
from docutils.parsers.rst import Directive


class BannerSmallDirective(Directive):
    """Python directive for sphinx banner_small"""
    required_arguments = 0
    option_spec = {
        'origin': str,
    }

    @staticmethod
    def html(origin: str) -> nodes.Node:
        # pylint: disable=line-too-long
        html = f'''
            <div class="cta-bar">
                <div class="cta-text">
                    ...
                </div>
                <button  onclick="window.open('https://console.inductiva.ai/api/register?utm_source=guide_{origin}', '_blank')" target="_blank" class="cta-button" id="login-btn-small">...</button>
            </div>
            '''
        return nodes.raw('', html, format='html')

    def run(self):
        origin = self.options.get('origin', '')
        return [BannerSmallDirective.html(origin)]
