# HSCM
**HSCM** (Higher-order Structural Centrality Methods) are centrality methods based on the higher-order structure of the network. This method extracts the higher-order structure of network, that is, the situation of nodes participating in the formation of cliques, and then cooperates with other imformation of the original network to weighting the base centrality methods. According to the different attributes considered, there are two concrete methods: CDR (Clustered Degree Rank with Higher-order Weighting) and CSR (Clustered Strength Rank with Higher-order Weighting).  
- **CDR**：Further highlights the importance of the node degree, and cooperates with the node clique structure information to weighting the base centrality score of the node.  
- **CSR**：Defines the strength of the node according to the situation of the node participating in the formation of cliques, and then cooperates with the clique structure information to weighting the base centrality score of the node.
## Dependency
This project runs in the environment of `python==3.10.12`.
## Installation
You can clone this project to your local machine in the following way:
`git clone -b master https://github.com/brleude/HOCR.git`
## Usage
We use an artificially constructed scale-free network for demonstration.The relevant data can be found in `HOCR/data/`.<br>**Run the following code in folder `HOCR/`  under the conda environment:**

```
conda create -n HOCR_test python=3.10.12
conda activate HOCR_test
pip install -r requirements.txt
python main.py \
-G data/artificial_network_edgelist.edgelist \
-b cc \
-o output

arguments:
-G, --network
		The path of network file in edgelist format.
-b, --base_scores
		Basic centrality scores, optional values are 'bc', 'cc',
		'dc', 'ec', 'pr' or a score dictionary.
-o, --output
		The path of output file.
```
**Other optional parameters：**
```
-k, --key_nodes
		List of key nodes or path, used for parameter tuning.
-n, --non_key_nodes
		List of non-key nodes or path, used for parameter tuning.
-m, --max_clique
		The size of the largest clique to be considered.
-r, --rank_type
		The method for parameter optimization, with options of
		cdr or csr, and the default value is csr.
-v, --verbose
		Whether to output the parameter optimization log, with
		options of True or False.

```
**result:** <br>You can find the result in the following file:
`output/output.txt`.<br>The output result is in the following format:
```
cdr: {'0': 0.005651006929750437, '1': 0.0020545118446591206, '2':0.0036751727019928787, '3': 0.01351130804370008,...}
csr: {'0': 0.011682150937732259, '1': 0.0020338003154715136, '2': 0.006527391739653339, '3': 0.057927416419532404,...}
```
The output result includes the CDR and CSR scores of all nodes in the network.

