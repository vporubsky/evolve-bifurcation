---


---

<h1 id="evolve-bifurcation">evolve-bifurcation</h1>
<h2 id="introduction">Introduction</h2>
<p>The evolve-bifurcation repository contains sourcecode for a  Pythonic bifurcation evolution algorithm, <a href="https://github.com/vporubsky/evolve-bifurcation/blob/master/evolveBifurcation.py">evolveBifurcation</a>. The module employs a global search algorithm, differential evolution, to iteratively search parameter space by mutating global parameters, floating species initial concentrations, and boundary species of a network model in antimony or sbml formats to optimize the system to contain eigenvalues characteristic of Hopf or turning point bifurcations.</p>
<h2 id="installation-instructions">Installation Instructions</h2>
<p>The evolveBifurcation module can be downloaded from <a href="https://github.com/vporubsky/evolve-bifurcation/blob/master/evolveBifurcation.py">evolveBifurcation.py</a></p>
<p>The user should download this python script to their working directory. After the download, the user can import the evolveBifurcation function:</p>
<pre><code>from evolveBifurcation import evolveBifurcation
</code></pre>
<p>After import, evolveBifurcation can now be called to evolve model parameters:</p>
<pre><code>optModel, fitness = evolveBifurcation(model='kholodenko.xml', ...
				 bifurType='oscillator')

</code></pre>
<h2 id="usage-examples">Usage Examples</h2>
<p>evolveBifurcation() can take model inputs that are an existing RoadRunner instance, an SBML model format (.xml), or an Antimony string model format. evolveBifurcation() returns the model as a RoadRunner instance that can be simulated with all parameters optimized for the assigned objective. The fitness value of the returned model is also returned. In the examples below, optModel receives the returned optimized model, and fitness receives the final objective function evaluation of the returned model. The model must always be specified by the user. All other input arguments have default values. The default bifurcation type to determine the search object is ‘oscillator’.</p>
<p>RoadRunner instance input example to optimize for a Hopf bifurcation (oscillator) with default settings of the tool:</p>
<pre><code>model = te.loadSBMLModel('kholodenko.xml')
optModel, fitness = evolveBifurcation(model, bifurType='oscillator')
</code></pre>
<p>SBML model (.xml) input example to optimize for a Hopf bifurcation (oscillator) with default settings of the tool:</p>
<pre><code>optModel, fitness = evolveBifurcation(model='kholodenko.xml', bifurType='oscillator')
</code></pre>
<p>Antimony string model (.ant) input example to optimize for a Hopf bifurcation (oscillator) with default settings of the tool:</p>
<pre><code>optModel, fitness = evolveBifurcation(model='relaxationOscillator.ant', bifurType='oscillator')
</code></pre>
<p>The returned optimized model can be simulated, as it is a functional RoadRunner instance:</p>
<pre><code>optModel, fitness = evolveBifurcation(model= 'kholodenko.xml', bifurType='oscillator', displayOpt = True)
simulationData = optModel.simulate(0, 1000, 1000)
optModel.plot() 
</code></pre>
<p>To search for a turning point bifurcations, change the bifurcation type input argument using the following command:</p>
<pre><code>bifurType = 'turningPoint'
</code></pre>
<p>The user can choose which parameters to optimize in their model by generating a list of the parameter names, entered as type str. The kholodenko.xml model contains 30 model parameters, floating species, and boundary species that would automatically be selected for optimization. To optimize only  J0_K1, J0_n, J7_KK8, and init([MKKK]), the user should add the following to the evolveBifurcation() input commands:</p>
<pre><code>parameterList = ['J0_K1', 'J0_n', 'J7_KK8', 'init([MKKK])']
</code></pre>
<p>Similarly, if the user wishes to generate custom parameter ranges that the algorithm will sample for each parameter, they may enter a list of tuple parameter ranges. The list must either have a length of 1, indicating that all parameters undergoing optimization will be subjected to this range of parameter values, or the list must have a length equal to the length of the parameterList. For example, if the user wishes to optimize the parameterList of 4 unique parameters specified above from the kholodenko.xml model described in the previous example, they should enter the following input command to evolveBifurcation():</p>
<pre><code>paramRanges = [(0, 10.0), (0.001, 1.0), (0.1, 5.0), (1.0, 100.0) ]
</code></pre>
<p>The maximum number of generations for differential evolution, the number of members in the population being optimized, the recombination constant value and the mutation constant value can all be altered to shift the algorithm from fast convergence to accurate convergence if desired. To alter these values, the user should add the following arguments to the list of input commands:</p>
<pre><code>maxGenerations = 500 
numMembers = 20 
mutationConst = 0.7 
recombinationConst = 0.7 
</code></pre>
<p>A threshold value will be used to determine whether the algorithm has converged upon a suitable solution. The threshold will either be based on the current fitness value, which is dependent on the evaluation of the bifurcation objective function, or on the eigenvalue of the member with the best fitness that best approximates the bifurcation behavior specified. The default is to use the fitness-based tolerance, with a value of 5, which has been empirically determined to be sufficient for most oscillatory models. These inputs can be altered by adding the following combinations to the input commands list:</p>
<pre><code>thresholdType = 'eigenvalue' 
threshold = 0.1 
</code></pre>
<p>While most problems work best without local minimization applied, the user can additionally implement a local minimization which uses a gradient-based method called the Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm. To do so, they should add the following to the list of input commands:</p>
<pre><code>localMin = True
</code></pre>
<p>The tool can show the progress of optimization of the differential evolution algorithm by progressively printing the generation and corresponding fitness value to the console by turning on the displayOpt input. To do so, add the following to the list of input commands:</p>
<pre><code>displayOpt = True
</code></pre>
<h2 id="authors">Authors</h2>
<ul>
<li><strong>Veronica Porubsky</strong>  -  <em>algorithm development, software development and programming, data production</em>  -  <a href="https://github.com/vporubsky">vporubsky</a></li>
<li><strong>Herbert Sauro</strong>  -  <em>algorithm development</em>  -  <a href="https://github.com/hsauro">hsauro</a></li>
</ul>
<p>See also the list of  <a href="https://github.com/vporubsky/evolve-bifurcation/graphs/contributors">contributors</a>  who participated in this project.</p>
<h2 id="license">License</h2>
<p>This project is licensed under the Apache License 2.0 - see the  <a href="https://github.com/vporubsky/evolve-bifurcation/blob/master/LICENSE">LICENSE.md</a>  file for details</p>
<h2 id="acknowledgments">Acknowledgments</h2>
<ul>
<li>Thank you to Kiri Choi, for valuable discussions on integrating steady state solvers and the functionality provided by Tellurium and libRoadRunner for implementing the bifurcation evolution algorithm.</li>
<li>Thank you to J. Kyle Medley for providing insight into the utility and shortcomings of several global and local optimization algorithms.</li>
</ul>

