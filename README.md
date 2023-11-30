# DrugRepurposing

Automated drug repurposing workflow for rare diseases using Flask.
A Masters project by Carmen Reep.

This drug repurposing app finds candidate compounds for a symptom of a specific rare disease.
For this, it uses a workflow that has four main steps:
1. Creation of an RDF knowledge graph by selecting a subnetwork from Monarch Initiative, adding drugs using the DGIdb database, creating drug-drug edges based on drug compound structure similarity.
2. Graph embedding using RDF2Vec.
3. Creation of edge representations (=training/test data) by fusing drug and gene feature vectors.
4. Training an XGBoost machine learning model which is evaluated using repeated stratified 10-fold cross validation, where the best model is used to predict unknown drug-gene interactions of interest.

The whole drug repurposing workflow is saved in five separate Python files: <code>Monarch.py</code>, <code>DGIdb.py</code>, <code>drugsimilarity.py</code>, <code>combine_graphs.py</code>, and <code>rdf2vec.py</code>. To run this whole process with Flask, three extra Python files were created: <code>\_\_init\_\_.py</code>, `routes.py`, and `forms.py`. 

The `__init__.py` file imports the `Flask` class from the `flask` package version 1.1.2 and instantiates the Flask app. The `routes.py` file specifies the routes of the app, which are acts that bind a URL to a view function (a function that responds to a request). The `forms.py` file includes form classes for user input and uses WTForms to render and validate the forms.

Each time the user makes a request (by clicking on a link or giving input via a form), the `render_template()` function is run, which generates output from a template file that is found in the \`templates\' folder of the app. For a nice general style of the app, each template file extends the `layout.html` template, which provides a site header and title and a section with links of all previously run predictions.

To start the server, the `run()` method of the `Flask` object is called. It returns the URL where the server is available. When opening the URL, the home page of the app is shown, which is rendered from the `home.html` template file. On this home page, there are several options for the user. The user can choose an existing disease graph, which has been previously created using a disease seed as input for Monarch, where the date of creation of the graph is specified. If the user is interested in a disease that is not in this list of previously created graphs, or the user wants a newer version of a graph, the user can create a new disease graph by specifying the phenotype MIM number of the disease of interest (e.g. \`143100\' for Huntington's disease). This enables the Python `Monarch.py` file to run the function that creates a new network from Monarch using this input number as seed. Creating a new network form Monarch takes some time, so the user is asked to be patient. After running a new disease graph, the graph is saved and added to the existing disease graphs on the homepage.

After selecting an existing disease graph, or creating a new disease graph, the `symptoms.html` template is rendered. Here, all symptoms of the disease of interest are listed, and the user is asked to select one symptom of interest. When a user selects a symptom, again, the Python `Monarch.py` file is called to run the function that creates a new network from Monarch, with this time only the chosen symptom as seed. After creating this symptom Monarch graph, this graph is merged with the disease Monarch graph to create the final Monarch graph. Then `DGIdb.py`, `drugsimilarity.py`, `combine_graphs.py`, and `rdf2vec.py` are run in this order. The last Python file saves the ranked predicted drug-gene interactions as an HTML file in the `templates` folder of the app. When everything is run, this predictions HTML file is rendered and the list of ranked predictions is shown on screen. This list includes the drugs, the genes each drug interacts with, together with the interaction type and confidence. Each drug and gene is presented with its label. Clicking on a label will direct the user to the URI link of the drug or gene, which shows all up-to-date information about that drug or gene.

On every page (home, symptom, prediction), we added the option to select a previously run prediction, to avoid long running times.

## Getting started
### Docker
Build docker image: <code>docker build --tag drugapp .</code>

Run the image: <code>docker run -d -p 5000:5000 --name drugapp drugapp</code>
### Without docker
Download the `app` folder, then run the Flask app with <code>python run.py</code>. 


