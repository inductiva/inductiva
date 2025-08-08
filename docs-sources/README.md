# Documentation Sources

This folder contains the **source files** for the documentation sections displayed on the [website](https://inductiva.ai/guides/documentation/cli/overview).

---

## New: Full Control Over the Navigation Menu

You no longer need to ask the dev team to add or change entries in the website's navigation menus.

All sections shown in the **"Simulators"** and **"Resources"** menus are now defined dynamically via a single JSON file:  
➡️ [`public/builds/menu.json`](../docs-sources/menu.json)

This gives you **full control** over:

- What appears in each section (e.g., Simulators, Resources)
- How many columns are used in the mega menu
- Whether to show or hide section titles

Just edit the `menu.json` structure, changes will be automatically reflected on the website, without any code changes needed.

---

## How It Works

Each subfolder in `docs-sources/` represents a distinct section of documentation, rendered on the website under the `/guides` path. For example:

- `docs-sources/openfast` → [`/guides/openfast`](https://inductiva.ai/guides/openfast)
- `docs-sources/how-it-works` → [`/guides/how-it-works`](https://inductiva.ai/guides/how-it-works)
- `docs-sources/cans` → [`/guides/cans`](https://inductiva.ai/guides/cans)
  
<br />

<img width="321" alt="Captura de ecrã 2025-03-26, às 17 35 13" src="https://github.com/user-attachments/assets/57abd61d-d7de-4bae-8760-73be342deb91" />
<img width="1500" alt="Captura de ecrã 2025-03-26, às 17 37 50" src="https://github.com/user-attachments/assets/f4829274-687a-4d03-af31-c41da1a8dcfd" />

You can view the documentation at these paths or through the **"Resources"** menu (if linked there).  

Previously, updating the Simulators or Resources menus required dev support. Now, everything is controlled via the menu.json file. **Now, it's fully configurable** through the `menu.json` file, offering total flexibility.

The **Markdown files** in this folder are copied to the website repository as part of the build process.  
On the website, **Sphinx** processes these Markdown files into HTML, which is then used to render the documentation.

A **[GitHub Action](https://github.com/inductiva/inductiva/blob/development/.github/workflows/website_deploy_trigger_on_docs_update.yml)** runs whenever changes are pushed to the `development` or `main` branches.  
This triggers a webhook that notifies the website repository to start a new deployment. During this process, the site pulls the latest Markdown, processes it with **Sphinx**, and updates the documentation.

---

## How to Edit the Documentation

To update documentation on the website, follow this process:

1. Make the necessary changes to the Markdown files inside `docs-sources/`.
2. To preview the changes on the **staging website**, create a pull request and then merge the changes to **`development`**.
3. To publish the changes to the **[production website](https://inductiva.ai/)**, create a pull request and then merge the changes to **`main`**.
4. If you want to **change what appears in the website's submenus** (like **Simulators** or **Resources**), edit the structure of the file `docs-sources/menu.json`.  
   This file defines which groups and items appear, how they are organized, and how many columns are shown.
