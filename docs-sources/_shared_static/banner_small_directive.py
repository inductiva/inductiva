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
        # pylint: disable=anomalous-backslash-in-string
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
                // Current URL query string, including '?'
                const params = new URL(window.parent.location.href).search;
                const parentPath = new URL(window.parent.location.href).pathname;
                
                // Replace "/" with "_", remove leading/trailing underscores if any
                const utmPath = parentPath.replace(/\//g, '_').replace(/^_+|_+$/g, '');

                const baseUrl = 'https://console.inductiva.ai/api/register?utm_cta_origin=guide_' 
                    + origin + '&utm_path=' + encodeURIComponent(utmPath);

                const url = params
                    ? baseUrl + '&' + params.slice(1)  // Remove the initial '?' and prepend '&'
                    : baseUrl;

                console.log("[Banner] Opening URL:", url);
                window.open(url, '_blank');
            }}
            </script>
            '''
        return nodes.raw('', html, format='html')

    def run(self):
        origin = self.options.get('origin', '')
        return [BannerSmallDirective.html(origin)]
