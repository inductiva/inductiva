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

            <div class="container">
        <div class="banner">
            <!-- Decorative star in top right -->
            <div class="star-icon">
                <svg width="32" height="32" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                    <path d="M9.937 15.5A2 2 0 0 0 8.5 14.063l-6.135-1.582a.5.5 0 0 1 0-.962L8.5 9.936A2 2 0 0 0 9.937 8.5l1.582-6.135a.5.5 0 0 1 .963 0L14.063 8.5A2 2 0 0 0 15.5 9.937l6.135 1.582a.5.5 0 0 1 0 .962L15.5 14.063a2 2 0 0 0-1.437 1.437l-1.582 6.135a.5.5 0 0 1-.963 0L9.937 15.5Z"/>
                </svg>
            </div>

            <!-- Main content -->
            <div class="content">
                <!-- Heading -->
                <h1 class="heading">Ready to dive in?</h1>

                <!-- Subtitle -->
                <p class="subtitle">
                    Get started for free today and earn 
                    <span class="highlight">$5 in credits</span>.
                </p>

                <!-- CTA Button -->
                <button class="cta-button">
                    Get Started
                    <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
                        <path d="M9.937 15.5A2 2 0 0 0 8.5 14.063l-6.135-1.582a.5.5 0 0 1 0-.962L8.5 9.936A2 2 0 0 0 9.937 8.5l1.582-6.135a.5.5 0 0 1 .963 0L14.063 8.5A2 2 0 0 0 15.5 9.937l6.135 1.582a.5.5 0 0 1 0 .962L15.5 14.063a2 2 0 0 0-1.437 1.437l-1.582 6.135a.5.5 0 0 1-.963 0L9.937 15.5Z"/>
                    </svg>
                </button>

                <!-- Feature points -->
                <div class="features">
                    <div class="feature-item">
                        <div class="feature-dot"></div>
                        <span class="feature-text">No credit card</span>
                    </div>
                    <div class="feature-item">
                        <div class="feature-dot"></div>
                        <span class="feature-text">Free forever</span>
                    </div>
                </div>
            </div>

            <!-- Decorative elements -->
            <div class="blur-element blur-bottom-left"></div>
            <div class="blur-element blur-top-right"></div>
        </div>
    </div>
            '''
        return [nodes.raw('', html, format='html')]
