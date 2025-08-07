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
                    <button onclick="openInductivaRegister('{origin}')" target="_blank" class="btn primary" id="login-btn-big" >
                    <span class="btn-main">...</span>
                    </button>
                </div>
                </div>
            </div>
            <script>
            function openInductivaRegister(origin) {{
                // current URL query string, including '?'
                const params = new URL(window.parent.location.href).search;
                const baseUrl = 'https://console.inductiva.ai/api/register?guides_cta_origin=guide_' + origin;
                
                const url = params
                    ? baseUrl + '&' + params.slice(1)  // Remove the initial '?' and prepend '&'
                    : baseUrl;

                console.log("[Banner] Opening URL:", url);
                window.open(url, '_blank');
            }}
            </script>
            <script>
            const button = document.getElementById('login-btn-big');

            function triggerShake() {{
                button.classList.add('shake');
                setTimeout(() => {{
                button.classList.remove('shake');
                }}, 500); // Match the animation duration
            }}

            setInterval(triggerShake, 10000);
            </script>
            '''
        return [nodes.raw('', html, format='html')]
