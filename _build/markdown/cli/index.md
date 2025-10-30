# Inductiva CLI

The **Inductiva CLI** helps you interact with the Inductiva API directly from the terminal, enabling you to manage computational resources, tasks, and much more.

This section documents every command and flag available in Inductiva’s command-line interface (CLI).

## CLI Command Overview

The table below shows the available CLI commands alongside their corresponding Python Client APIs and documentation guide. Most functionality is available through both interfaces, allowing you to choose the tool that best fits your workflow. The CLI is ideal for quick operations, while the Python Client offers a more programmatic control and integration with your code.

For detailed guidance on when to use each interface and how they work together, see our [Interfaces with the API](http://inductiva.ai/guides/how-it-works/building-blocks/index) guide.

| Command                         | Python Client                                                                        | Resource Guide                                                                             |
|---------------------------------|--------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
| [`auth`](auth.md)               | –                                                                                    | [Authentication Guide](https://inductiva.ai/guides/how-it-works/get-started/install-guide) |
| [`user`](user.md)               | [Users API](https://inductiva.ai/guides/api-functions/api/inductiva.users)           | –                                                                                          |
| [`tasks`](tasks.md)             | [Tasks API](https://inductiva.ai/guides/api-functions/api/inductiva.tasks)           | [Tasks Guide](https://inductiva.ai/guides/how-it-works/tasks/index)                        |
| [`task-runner`](task-runner.md) | –                                                                                    | [BYOH Guide](https://inductiva.ai/guides/expand/use-local-task-runner/index)               |
| [`projects`](projects.md)       | [Projects API](https://inductiva.ai/guides/api-functions/api/inductiva.projects)     | [Projects Guide](https://inductiva.ai/guides/scale-up/projects/index)                      |
| [`storage`](storage.md)         | [Storage API](https://inductiva.ai/guides/api-functions/api/inductiva.storage)       | [Storage Guide](https://inductiva.ai/guides/how-it-works/intro/data_flow)                  |
| [`resources`](resources.md)     | [Resources API](https://inductiva.ai/guides/api-functions/api/inductiva.resources)   | [Resouces Guide](https://inductiva.ai/guides/how-it-works/machines/index)                  |
| [`simulators`](simulators.md)   | [Simulators API](https://inductiva.ai/guides/api-functions/api/inductiva.simulators) | [Simulators Guide]()                                                                       |
| [`containers`](containers.md)   | –                                                                                    | [BYOS Guide](https://inductiva.ai/guides/expand/bring-your-own-software/index)             |

### Notes

- Use `--help` with any command for more detailed usage and available flags.
- For safety, commands like `resources terminate` and `containers remove` include a `--yes` flag to bypass confirmation prompts.

---

## Set Up & Authentication

The CLI is automatically installed when you install
Inductiva’s [Python Client](../api/index.md). However, before using the CLI,
you need to authenticate with your Inductiva account.

```sh
inductiva auth login
```

When prompted, copy paste your personal API key that is availalbe
from your [User Account](https://console.inductiva.ai/account/profile)
page in the Web Console.

```sh
     ___  _   _  ____   _   _   ____  _____  ___ __     __ _
    |_ _|| \ | ||  _ \ | | | | / ___||_   _||_ _|\ \   / // \
     | | |  \| || | | || | | || |      | |   | |  \ \ / // _ \
     | | | |\  || |_| || |_| || |___   | |   | |   \ V // ___ \
    |___||_| \_||____/  \___/  \____|  |_|  |___|   \_//_/   \_\
    
    To log in, you need an API Key. You can obtain it from your account at https://console.inductiva.ai/account.
    Please paste your API Key here: 
```

## Need Help?

Use the `--help` or `-h` flag with **any command** for more details:

```sh
inductiva tasks --help
```


            <div class="banner">
                <div class="banner-content">
                <div class="text">
                    <p class="headline">...</p>
                    <p class="subtext">...</p>
                </div>
                <div class="buttons">
                    <button onclick="openInductivaRegister('cli')" target="_blank" class="btn primary" id="login-btn-big" >
                    <span class="btn-main">...</span>
                    </button>
                </div>
                </div>
            </div>
            <script>
            function openInductivaRegister(origin) {
                // Current URL query string, including '?'
                const params = new URL(window.parent.location.href).search;
                const parentPath = new URL(window.parent.location.href).pathname;
                
                // Replace "/" with "_", remove leading/trailing underscores if any
                const utmPath = parentPath.replace(/\\//g, '_').replace(/^_+|_+$/g, '');

                // Get referrer domain only, remove protocol and www
                let referrerDomain = '';
                try {
                    const refUrl = new URL(window.parent.document.referrer);
                    referrerDomain = refUrl.hostname.replace(/^www\\./, ''); // Remove www.
                } catch (e) {
                    // If referrer is empty or invalid, leave as empty string
                    referrerDomain = '';
                }

                // Sanitize: allow only alphanumerics, "-", "_", "."
                const utmReferrer = encodeURIComponent(referrerDomain.replace(/[^a-zA-Z0-9_\\-\\.]/g, '_'));

                const baseUrl = 'https://console.inductiva.ai/api/register?utm_cta_origin=guide_' 
                    + origin + '&utm_path=' + encodeURIComponent(utmPath)
                    + '&utm_ref=' + utmReferrer;

                const url = params
                    ? baseUrl + '&' + params.slice(1)  // Remove the initial '?' and prepend '&'
                    : baseUrl;

                console.log("[Banner] Opening URL:", url);
                window.open(url, '_blank');
            }
            </script>
            <script>
            const button = document.getElementById('login-btn-big');

            function triggerShake() {
                button.classList.add('shake');
                setTimeout(() => {
                button.classList.remove('shake');
                }, 500); // Match the animation duration
            }

            setInterval(triggerShake, 10000);
            </script>
