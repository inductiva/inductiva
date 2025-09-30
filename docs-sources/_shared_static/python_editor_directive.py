"""Python directive for sphinx call-to-action button"""
from docutils import nodes
from docutils.parsers.rst import Directive


class PythonEditorDirective(Directive):
    """Python directive for sphinx call-to-action button"""
    required_arguments = 0
    option_spec = {
        'origin': str,
        'text': str,
        'subtext': str,
        'button_text': str,
        'url': str,
    }

    def run(self):
        origin = self.options.get('origin', '')
        text = self.options.get(
            'text',
            'Try our online Python Editor – run simulations in your browser')
        button_text = self.options.get('button_text', 'Start Simulating Now')
        base_url = self.options.get('url',
                                    'https://console.inductiva.ai/editor')

        # pylint: disable=line-too-long
        # pylint: disable=anomalous-backslash-in-string
        html = f'''
            <div style="text-align: center; margin: 20px 0;">
                <div style="font-size: 12px; margin-bottom: 6px;">{text}</div>
                <a href="#" 
                   onclick="openInductivaEditor('{base_url}', '{origin}'); return false;"
                   style="display: inline-block; width: 55%; padding: 16px 24px; font-size: 14px; font-weight: bold; background-color: var(--playground-button); color: black; text-decoration: none; text-align: center; border-radius: 8px;">
                    {button_text}
                </a>
            </div>
            <script>
            function openInductivaEditor(baseUrl, origin) {{
                // Current URL query string, including '?'
                const params = new URL(window.parent.location.href).search;
                const parentPath = new URL(window.parent.location.href).pathname;
                
                // Replace "/" with "_", remove leading/trailing underscores if any
                const utmPath = parentPath.replace(/\//g, '_').replace(/^_+|_+$/g, '');

                // Get referrer domain only, remove protocol and www
                let referrerDomain = '';
                try {{
                    const refUrl = new URL(window.parent.document.referrer);
                    referrerDomain = refUrl.hostname.replace(/^www\./, ''); // Remove www.
                }} catch (e) {{
                    // If referrer is empty or invalid, leave as empty string
                    referrerDomain = '';
                }}

                // Sanitize: allow only alphanumerics, "-", "_", "."
                const utmReferrer = encodeURIComponent(referrerDomain.replace(/[^a-zA-Z0-9_\-\.]/g, '_'));

                // Build URL with UTM parameters
                const url = new URL(baseUrl);
                
                // Add utm_cta_origin parameter
                if (origin) {{
                    url.searchParams.set('utm_cta_origin', 'guide_' + origin);
                }}
                
                // Add utm_path parameter
                if (utmPath) {{
                    url.searchParams.set('utm_path', utmPath);
                }}
                
                // Add utm_ref parameter
                if (utmReferrer) {{
                    url.searchParams.set('utm_ref', utmReferrer);
                }}
                
                // Preserve existing URL parameters
                if (params) {{
                    const existingParams = new URLSearchParams(params);
                    existingParams.forEach((value, key) => {{
                        // Don't override our UTM parameters
                        if (!url.searchParams.has(key)) {{
                            url.searchParams.set(key, value);
                        }}
                    }});
                }}

                console.log("[Call to Action] Opening URL:", url.toString());
                window.open(url.toString(), '_blank');
            }}
            </script>
            '''
        return [nodes.raw('', html, format='html')]
