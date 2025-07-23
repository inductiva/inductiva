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
                    <span class="btn-main">...</span>
                    </button>
                </div>
                </div>
            </div>
            <script>
            const button = document.getElementById('login-btn-big');

            function triggerShake() {{
                button.classList.add('shake');
                setTimeout(() => {{
                button.classList.remove('shake');
                }}, 500); // Match the animation duration
            }}

            // Shake every 5 seconds
            setInterval(triggerShake, 5000);
            </script>
            '''
        return [nodes.raw('', html, format='html')]
