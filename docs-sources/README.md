<!-- # Documentation Deployment Workflow

This directory contains the source files for the documentation presented on the website. -->
<!-- To keep the documentation up to date in the **staging environment**, we use a [GitHub Action workflow](https://github.com/inductiva/website-new/actions/workflows/deploy-docs-to-staging.yml) that automatically deploys changes to the `development` branch. -->

<!-- ### How to Trigger the Deployment To Development/Staging

To deploy changes from any branch to the `development` environment, follow these steps:

1. **Make Changes**:  
   Make any necessary changes in the `docs-sources` folder. This folder contains the source files for the documentation.

2. **Trigger the Workflow**:  
   - Go to the **GitHub repository** page.
   - Navigate to the **Actions** tab.
   - Find and select the workflow named **Deploy Docs to Staging (via development branch)**.
   - Under the **"Run workflow"** dropdown, choose the branch containing the changes you want to deploy (e.g., docs-update).
   - Click **"Run workflow"** to trigger the process.
  
     <img width="1490" alt="Captura de ecrã 2025-03-26, às 17 05 43" src="https://github.com/user-attachments/assets/df2b611e-437a-44cc-ba84-5f7eac067303" />

3. **What Happens Next**:  
   - The workflow will check if there are any changes in the `docs-sources` folder.
   - If changes are found, it will automatically **commit** them to the `development` branch and **push** the changes to the repository.
   - This push to `develipment` will trigger the **Google Cloud Build** process, deploying the updated documentation to the development environment (this can take a while — a few minutes).

### Why Use This Workflow?

This workflow allows automatic updates to the documentation without needing to go through pull requests or manual deployments. It ensures that only changes to the documentation (in the `docs-sources` folder) are deployed, keeping the rest of the code untouched.

### Important Notes

- **Changes outside the `docs-sources` folder will not trigger the workflow**. If you modify files outside of this folder, the workflow will not take any action, and those changes will not be deployed.
  
- If you need to **add or remove documentation options** presented the submenu under the **Resources** section that shows some of these docs, this change must be made **via code**. For such modifications, please **contact Sara** to handle those changes. -->

### How to View the Documentation Changes

As you make changes to the documentation, they are reflected in the website, depending on the structure of the `docs-sources` folder. Each subfolder within `docs-sources` represents a distinct root of the documentation, and these sections can be accessed under the **/guides** path on the website.

For example, let's use [**openfast**](https://website-staging.inductiva.ai/guides/openfast) documentation as an example:

- **`docs-sources/openfast` → `/guides/openfast`**
- `docs-sources/tutorials` → `/guides/tutorials`
- `docs-sources/reef3d` → `/guides/reef3d`

<img width="321" alt="Captura de ecrã 2025-03-26, às 17 35 13" src="https://github.com/user-attachments/assets/2ba1c174-19c6-4a0c-a49d-bf66ecd1e08f" />

<img width="1500" alt="Captura de ecrã 2025-03-26, às 17 37 50" src="https://github.com/user-attachments/assets/9e1d7b8a-19ec-43f2-af82-03e179083f04" />

### What Happens When You Make Changes

The documentation files within the `docs-sources` folder are written in **Markdown**. When changes are made, **Sphinx** is used to process these Markdown files and generate the corresponding HTML files. These HTML files are then rendered by the website.

### What to Do After Deploying

Once the deployment is triggered and the changes are pushed to the `development` branch, you can view the updated pages by navigating to the **[staging environment](https://website-staging.inductiva.ai/)** and going to the `/guides/${doc-folder-name}`. Each modified folder in `docs-sources` will correspond to a specific path under **/guides**.

This way, you can easily preview the changes to each section as they are deployed to the staging environment.

<!-- ---

## For Developers: How to View Changes Locally

If you want to test the changes locally before deploying them to the staging environment, follow these steps:

### 1. **Install Dependencies**

```sh
pip install -r docs-sources/requirements.txt  # Install Sphinx dependencies
npm i  # Install project dependencies
```

### 2. **Start the Project Locally**

```sh
npm run dev  # Start the project in development mode
```

If you are still having trouble initiating the project, check [here](https://github.com/inductiva/website-new/blob/development/README.md) to ensure that you have set up everything necessary for the Nuxt project.

### 3. **After Making Changes to the Docs**

Whenever you make changes to the documentation files and want to see them reflected locally, run:

```sh
npm i
npm run dev
```

or

```sh
npm i && npm run dev
```

This will build the new changes and restart the project to ensure that updates are visible. -->
