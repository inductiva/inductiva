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
                    <button onclick="handleClick()" target="_blank" class="btn primary" id="login-btn-big" >
                    <span class="btn-main">Get Started</span>
                    </button>
                </div>
                </div>
            </div>

            <script>
            function handleClick() {{
            // Give GA4 and GTM time to fire
            setTimeout(function() {{
                window.open('https://console.inductiva.ai/?utm_source=guide_{origin}&utm_medium=button&utm_campaign=signup', '_blank');
            }}, 500); // Timeout for google tags to do some work
            }}
            </script>
            '''
        return [nodes.raw('', html, format='html')]
