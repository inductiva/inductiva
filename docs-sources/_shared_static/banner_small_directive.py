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
                <div class="cta-text" style="text-align: center">
                    ...
                </div>
                <button
                    onclick="openInductivaRegister('{origin}')"
                    class="cta-button"
                    id="login-btn-small">...</button>
            </div>

            <script>
            function openInductivaRegister(origin) {{
                // current URL query string, including '?'
                const params = new URL(window.parent.location.href).search;
                const baseUrl = 'https://console.inductiva.ai/api/register?guides_cta_origin=guide_' + origin;
                
                const url = params
                    ? baseUrl + '&' + params.slice(1)  // Remove the initial '?' and prepend '&'
                    : baseUrl;

                console.log("[BannerSmall] Opening URL:", url);
                window.open(url, '_blank');
            }}
            </script>
            '''
        return nodes.raw('', html, format='html')

    def run(self):
        origin = self.options.get('origin', '')
        return [BannerSmallDirective.html(origin)]
