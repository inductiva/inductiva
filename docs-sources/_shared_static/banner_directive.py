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
        # pylint: disable=anomalous-backslash-in-string
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
