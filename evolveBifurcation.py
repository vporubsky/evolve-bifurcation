import tellurium as te
import tesbml
from roadrunner import RoadRunner
import numpy as np
import sys
def evolveBifurcation(model, bifurType = 'oscillator', parameterList = None, paramRanges = None, maxGenerations = 100, numMembers = 50, mutationConst = 0.5, recombinationConst = 0.5, thresholdType = 'fitness', threshold = 0.1, localMin = False, displayOpt = False):    
    """
    Author: Veronica Porubsky

    evolveBifurcation(model, bifurType = 'oscillator', parameterList = None, paramRanges = None, maxGenerations = 100, numMembers = 50, mutationConst = 0.5, recombinationConst = 0.5, thresholdType = 'fitness', threshold = 0.1, localMin = False, displayOpt = False)
    ======================
    
    The module will evolve global parameters, floating species initial 
    concentrations, and boundary species of the loaded antimony or sbml model
    to optimize for eigenvalues characteristic of Hopf or turning point 
    bifurcations.
    
    Inputs:
        
    model : executable model string, SBML (.xml) or Antimony (.ant) model 
        string file evolveBifurcation takes models in the form of Antimony or 
        SBML strings which have either an existing instance or which are saved 
        independently as files with the extension .ant or .xml and entered as 
        type str
             
    bifurType : str
        The bifurcation type that the user wishes to optimize for - either for 
        Hopf bifurcations to evolve parameters which produce an oscillatory
        output or for a turning point bifurcation which will attempt to evolve 
        parameters that may generate bistability. Should select either:
            - 'oscillator'
            - 'turningPoint'

    parameterList : list
        parameterList specifies the parameters to be optimized. If the model 
        contains parameters named k1, k3_deg, and k5, which the user wishes to
        optimize, these should be entered as follows: 
        parameterList = ['k1', 'k3_deg', 'k5']. The default value is 
        parameterList = None, in which the algorithm automatically optimizes
        all global parameters and all floating and boundary species initial
        conditions.
        
    paramRanges : sequence
        paramRanges imposes a minimum and maximum on the range of values 
        that the parameters being optimized can be assigned. Should be entered 
        as follows: paramRanges = (lower bound, upper bound). The lower and 
        upper bound should both of type float.
        
    maxGenerations : int
        maxGenerations sets the determines the maximum number of generations 
        the differential evolution algorithm will complete.
    
    numMembers : int
        Specifies the number of members included in the population for 
        differential evolution. Increasing the number of members may slow
        convergence but improve the quality of the solution.
    
    mutationConst : float
        The mutationConst parameter specifies the scaling factor during 
        
    recombinationConst : float
        The recombinationConst is used to determine whether or not 
        a parameter value in the mutant trial vector will undergo mutation or 
        whether it will accept the value of the original member.
    
    thresholdType : str
        The thresholdType determines whether the specified threshold should be
        checked with respect to the best fitness value of the differential
        evolution algorithm or to the smallest eigenvalue real component 
        magnitude. Should select either:
            - 'fitness'
            - 'eigenvalue'
    
    threshold : float
        The threshold is used to fully define stopping criteria for
        the differential evolution algorithm.
    
    localMin : bool
        Setting this argument to True will allow the optimization algorithm to 
        perform a final local minimization using the 
        Broyden-Fletcher-Goldfarb-Shanno algorithm, maintaining the bounded 
        parameter ranges used during global minimization with differential 
        evolution.
    
    displayOpt : bool
        Use diplayOpt = True to display the progress of optimization. Will 
        print to the console the optimized objective function evaluation at the 
        current generation of differential evolution or iteration of BFGS-B. 
    
    Returns:
        
    bifurcationModel : executable model/ RoadRunner object
        The resulting model with all differential evolution-optimized 
        parameters to produce the specified bifurcation type. This can be 
        accessed from the first element of the evolveBifurcation return, and 
        can be used directly for simulation in tellurium. 
        
    fitness : float
        Returns the minimized objective function evaluation, which is the 
        fitness value for the returned model.
    
    :Example:
        
    >>> import tellurium
    >>> from bifurcationEvolution import evolveBifurcation
    >>> r, funEval = evolveBifurcation('oscillatoryModel.xml', 'oscillator') # assign model with evolved parameters to roadrunner object 'r'
    >>> r.simulate() # simulate the updated model which is now available as a roadrunner object 
    """
    try:
        if type(model) == str:
            if model.endswith('.xml'):
                try:
                    bifurcationModel = te.loadSBMLModel(model)
                    bmSBMLPP = bifurcationModel.getParamPromotedSBML(bifurcationModel)
                    bifurcationModel = te.loadSBMLModel(bmSBMLPP)
                except:
                    print('Could not promote local parameters to global parameters. Continuing optimization without parameter promotion.')
                    bifurcationModel = te.loadSBMLModel(model)
            elif model.endswith('.ant'):
                try:
                    bifurcationModel = te.loada(model)
                    bmSBML = bifurcationModel.getCurrentSBML()
                    bmSBMLPP = bifurcationModel.getParamPromotedSBML(bmSBML)
                    bifurcationModel = te.loadSBMLModel(bmSBMLPP)
                except:
                    print('Could not promote local parameters to global parameters. Continuing optimization without parameter promotion.')
                    bifurcationModel = te.loada(model)
        elif isinstance(model, RoadRunner):
            try:
                bifurcationModel = model
                bmSBML = bifurcationModel.getCurrentSBML()
                bmSBMLPP = bifurcationModel.getParamPromotedSBML(bmSBML)
                bifurcationModel = te.loadSBMLModel(bmSBMLPP)
            except:
                print('Could not promote local parameters to global parameters. Continuing optimization without parameter promotion.')
                bifurcationModel = model
        else:
            raise RuntimeError('Input not supported: pass an exisiting roadrunner.RoadRunner instance, Antimony (.ant) file, or an SBML (.xml) file.')
    except Exception as e: print(e)
    if not bifurType in ['oscillator', 'turningPoint']:
        raise RuntimeError("bifurType must be either 'oscillator' or 'turningPoint'.")
    bifurcationModel.reset()
    bifurcationModel.conservedMoietyAnalysis = True
    parameterList, paramValues  = generateParameterList(bifurcationModel, parameterList)
    if paramRanges:
        if len(paramRanges) == 1:
            paramRanges = paramRanges*len(parameterList)
    else:
        paramRanges = generateParameterRanges(paramValues)
    toBifurcationObjFunc = (bifurcationModel, bifurType, parameterList)
    try:
        bifurcationOptParams, fitness = differentialEvolution(bifurcationObjFunc, parameterRanges = paramRanges, maxGenerations = maxGenerations, displayDEProgress = displayOpt, populationSize = numMembers, crossoverProbability = recombinationConst, mutationConstant = mutationConst, thresholdType = thresholdType, threshold = threshold, bfgsMin = localMin, differentialEvArguments = toBifurcationObjFunc) 
    except Exception as e: print(e)
    bifurcationModel = setModelParams(bifurcationModel, parameterList,  bifurcationOptParams)
    return bifurcationModel, fitness

def generateParameterList(bifurcationModel, parameterList=None):
    """
    Generates a list of global parameters, floating and boundary species, and 
    a second list containing the current values for these parameters in the 
    loaded model. Conserved sum parameters, generated becaused 
    conserved moiety analysis is employed, are removed. Parameters defined
    by an assignement rule are removed.
    """
    if not parameterList:
        parameterList = bifurcationModel.getGlobalParameterIds() + list(bifurcationModel.getFloatingSpeciesInitialConcentrationIds()) + bifurcationModel.getBoundarySpeciesIds()          
        paramValues = list(bifurcationModel.getGlobalParameterValues()) + list(bifurcationModel.getFloatingSpeciesConcentrations()) + list(bifurcationModel.getBoundarySpeciesConcentrations())
    else:
        paramValues = []
        for parameter in parameterList:
            paramValues.append(bifurcationModel.getValue(parameter))
    doc = tesbml.readSBMLFromString (bifurcationModel.getCurrentSBML())
    model = doc.getModel();
    lenListRules = len(model.getListOfRules())
    assignmentRuleParams = []
    for i in range(lenListRules):
        assignmentRule = model.getRule(i)
        if (assignmentRule.getTypeCode() == tesbml.SBML_ASSIGNMENT_RULE):
            assignmentRuleParams.append(assignmentRule.getVariable())
    for element in range(len(parameterList)):
        for n in assignmentRuleParams:
            if (parameterList[element].startswith('_CSUM')) or (parameterList[element] == n):
                parameterList[element] = 'remove'
                paramValues[element] = 'remove'
    parameterList = list(filter(lambda a: a != 'remove', parameterList))
    paramValues = list(filter(lambda a: a != 'remove', paramValues))
    return parameterList, paramValues

def generateParameterRanges(paramValues):
    """
    Generates a sequence containing bounded ranges for each parameter using the 
    current parameter values from the model passed to evolveBifurcation. Only 
    called if parameter ranges are not specified by the user. 
    """
    paramRanges = [()]*len(paramValues)
    for i in range(len(paramValues)):
        if paramValues[i] == 0.0:
            paramRanges[i] = (1E-25, 10.0)
        elif paramValues[i] <= 10.0 and not paramValues[i] < 0.0:
            paramRanges[i] = (paramValues[i]/10.0, 10.0)
        elif paramValues[i] > 10.0:
            paramRanges[i] = (paramValues[i]/10.0, 2*paramValues[i])
    return paramRanges

def bifurcationObjFunc(parameterValues, *toBifurcationObjFunc):
    """
    Returns the fitness value of the parameter values passed to the function, 
    corresponding to optimization for Hopf and turning point bifurcations. The 
    fitness is dependent on the eigenvalues of the loaded model. If the loaded 
    model cannot achieve steady state with the parameter values (s) passed to 
    the function, or if the routine raises an error, the fitness is assigned a 
    sufficiently large value (1E20), which ensures that the member is excluded 
    from the population undergoing optimization within the differential evolution.
    """
    try:
        bifurcationModel, bifurType, parameterList = toBifurcationObjFunc[0:3]
        bifurcationModel = setModelParams(bifurcationModel, parameterList, parameterValues)
        try:
            realComponent, imaginaryComponent = getSSEigenvalComponents(bifurcationModel)
        except:
            fitness = 1E20
            return fitness
        realFit = 1.0 
        imagFit = 1.0
        penalty = 0.0
        if bifurType == 'oscillator': 
            for n in range(len(realComponent)):
                if imaginaryComponent[n]:
                    realFit = np.abs(realFit*realComponent[n])
                imagFit = imagFit*(1.0 - ((0.99)*np.exp(-1.0*np.abs(imaginaryComponent[n])))) 
                if np.abs(realComponent[n]) < 10E-6 and np.abs(imaginaryComponent[n]) < 10E-3: 
                    penalty = penalty + 1E20
            fitness = (realFit/imagFit) + penalty
            return fitness
        elif bifurType == 'turningPoint':
            hasRealEigenValue = False
            bufferArray = []
            realFit = realComponent[0]
            imaginaryFit = imaginaryComponent[0]
            for i in range(1, len(realComponent)):
                tempRealFit = realFit
                realFit = realFit*realComponent[i] - imaginaryFit*imaginaryComponent[i]
                imaginaryFit = tempRealFit*imaginaryComponent[i] + imaginaryFit*realComponent[i]
            eigenValsProd = np.abs(realFit)
            for i in range(len(realComponent)):
                if np.abs(imaginaryComponent[i]) > 0.0:
                    bufferArray = np.append(bufferArray, 1E6)
                else:
                    bufferArray = np.append(bufferArray, np.abs(realComponent[i]))
                    hasRealEigenValue = True
            if hasRealEigenValue == True:
                minRealEigenValueIdx = np.argmin(bufferArray) 
                realFit = 1.0
                imaginaryFit = 0.0
                for i in range(len(realComponent)):
                    if i != minRealEigenValueIdx:
                        realFit = realFit*realComponent[i] - imaginaryFit*imaginaryComponent[i]
                        imaginaryFit = tempRealFit*imaginaryComponent[i] + imaginaryFit*realComponent[i]                        
                if len(realComponent) == 1 and realComponent[0] != 0.0: 
                    realEigenValsExceptMinProd = 0.0
                else:
                    realEigenValsExceptMinProd = np.abs(realFit)
            else:
                realEigenValsExceptMinProd = eigenValsProd
            denomProd = (1.0 - ((0.99)*np.exp(-1.0*realEigenValsExceptMinProd)))     
            fitness = eigenValsProd/denomProd
            return fitness
    except:
        fitness = 1E20
        return fitness

def getSSEigenvalComponents(model):
    """
    Returns the real and complex components of model eigenvalues at steady state.
    """
    model = steadyState(model)
    eigenVals = model.getReducedEigenValues()
    return np.real(eigenVals), np.imag(eigenVals)

def setModelParams(model, parameterList, parameterValues):
    """
    Sets the model parameters to the values passed to the function.
    """
    for idx, parameter in enumerate(parameterList):
        model.setValue(parameter, parameterValues[idx])
    return model

def steadyState(model, tolerance = 1E-15, maxIterations = 500):
    """
    Brings the model to steady state.
    """
    for i in range(maxIterations): 
        jacobianInverse = np.linalg.inv(model.getReducedJacobian())
        direction = -np.dot(jacobianInverse , model.dv()[0])
        stepSizeMultiplier = 1.0
        for idx, floatingSpecies in enumerate(model.getIndependentFloatingSpeciesIds()):
            numIter = 0
            model[floatingSpecies] = model[floatingSpecies] + stepSizeMultiplier*direction[idx]
            while model[floatingSpecies] < 0.0:
                stepSizeMultiplier = stepSizeMultiplier/2
                model[floatingSpecies] = model[floatingSpecies] + stepSizeMultiplier*direction[idx]
                if numIter > 250:
                    raise RuntimeError
                numIter += 1
        if np.linalg.norm(model.dv()) < tolerance:
            return model
    raise RuntimeError

def differentialEvolution(objFunc, parameterRanges, maxGenerations = 100, crossoverProbability = 0.5, mutationConstant = 0.5, populationSize = 40, thresholdType = 'fitness', threshold = 0.1, bfgsMin = False, displayDEProgress = False, differentialEvArguments = ()):
    """
    Generates a population of members that represent distinct locations in 
    parameter space, where, if there are N parameters, parameter space contains 
    N dimesions. Each member of the population then undergoes recombination, 
    mutation, and selection opeerations as the best fitness is continuosly 
    updated while the algorithm approaches a suitable solution. 
    
    Recombination involves selecting a random number and comparing it to the 
    recombination constant - if the number is smaller than the recombination
    constant, the current parameter value in the list undergoes mutation. If the 
    random number is larger than the recombination constant, the parameter value
    in the list is obtained from the current population member.
    
    Mutation involves selecting three members of the population which won a 
    single round of tournament selection. Members being compared in a round
    of tournament selection must be unique, however, the three 
    members involved in mutation may correspond to the same member. Mutation
    is an arithmetic operation where the mutated element is a sum of a single 
    best mutation member element and the difference between the remaining two
    mutation member elements multiplied by a scaling factor.
    
    Selection involves comparing the fitness values of the original population 
    member and mutated member. The member with a lower fitness is passed back
    into the population. 
    
    Termination of the algorithm is reached when the maximum number of 
    generations is exceeded, when the fitness value is not changing significantly
    over several generations, or when the threshold is reached. The user can 
    specify the type of theshold used through evolveBifurcation - 'fitness' simply 
    requires that the fitness value drops below a defined threshold, 'eigenvalue'
    requires that the smallest real component eigenvalue magnitude drops below
    a defined threshold.
    """
    bifurcationModel, bifurType, parameterList = differentialEvArguments[0:3]
    populationArray = [[] for member in range(populationSize)]
    populationFitness = np.zeros(populationSize)
    rejectPop = 0
    for i in range(populationSize):
        acceptMember = False
        while acceptMember == False:
            newMember = np.zeros(len(parameterRanges))
            for j in range(len(parameterRanges)):
                if parameterRanges[j][1] == 10.0:
                    newMember[j] = np.exp(np.random.uniform(np.log(parameterRanges[j][0]), np.log(parameterRanges[j][1])))
                else:
                    newMember[j] = np.random.uniform(parameterRanges[j][0], parameterRanges[j][1])
            populationFitness[i] = objFunc(newMember, *differentialEvArguments)
            if populationFitness[i] >= 1E20 or np.isnan(populationFitness[i]): 
                if rejectPop > 50000:
                    sys.exit('Unable to construct starting population of accepted members.')
                rejectPop += 1
                continue
            else:
                populationArray[i] = newMember 
                acceptMember = True
    numGenerations = 0
    terminationReached = False
    while terminationReached == False:
        for currentIdx in range(populationSize):
            mutationAgent1, mutationAgent2, mutationAgent3 = generateMutationAgents(populationArray, tournamentSelection(populationFitness))
            mutatedMember = np.zeros(len(parameterRanges))
            for i in range(len(parameterRanges)):
                randomNumber = np.random.uniform(0.0, 1.0)
                if randomNumber < crossoverProbability:
                    mutatedMember[i] = mutationAgent1[i] + mutationConstant*(mutationAgent2[i] - mutationAgent3[i]) 
                    if parameterRanges[i][0] <= mutatedMember[i] and mutatedMember[i] <= parameterRanges[i][1]:
                        continue
                    else:
                        mutatedMember[i] = populationArray[currentIdx][i]        
                else:
                    mutatedMember[i] = populationArray[currentIdx][i]
            currentMemberFitness = populationFitness[currentIdx]
            mutatedFitness = objFunc(mutatedMember, *differentialEvArguments)
            if not currentMemberFitness < mutatedFitness:
                populationArray[currentIdx] = mutatedMember
                populationFitness[currentIdx] = mutatedFitness
        bestMemberIdx = np.argmin(populationFitness)
        bestFitness = populationFitness[bestMemberIdx]
        bestModel = setModelParams(bifurcationModel, parameterList, populationArray[bestMemberIdx])
        realComponent, imaginaryComponent = getSSEigenvalComponents(bestModel)
        numGenerations += 1
        if displayDEProgress == True:
            print('Differential evolution generation ', str(numGenerations), ' fitness: f(x) = ', str(bestFitness))
        terminationReached = checkDETerminationCondition(numGenerations, maxGenerations, thresholdType, threshold, bestFitness, populationFitness, realComponent, imaginaryComponent, bifurType)
    if bfgsMin:
        return BFGS(objFunc, initialGuess = populationArray[bestMemberIdx], parameterRanges = parameterRanges, tol = 1E-25, displayBFGSProgress = displayDEProgress, BFGSArguments = differentialEvArguments) 
    return populationArray[bestMemberIdx], bestFitness

def generateMutationAgents(populationArray, agentIndices):
    """
    Returns three members from the population which won a single round of 
    tournament selection for mutation. 
    """
    bestAgents={}
    for i in range(0, 3):
        bestAgents["best{0}".format(i)] = populationArray[int(agentIndices[i])]
    return bestAgents["best0"], bestAgents["best1"], bestAgents["best2"]    

def tournamentSelection(populationFitness):
    """
    Compares the fitness values which correspond to two distinct population members 
    from the differential evolution population and returns the index for
    the position of the member with better fitness. Three indices are returned,
    corresponding to three winning members of a single round of tournament 
    selection. 
    """
    agentIndices=[]
    for i in range(1, 4):
        idx1, idx2 = chooseRandomMembers(len(populationFitness))
        if populationFitness[idx1] <= populationFitness[idx2]:
            agentIndices = np.append(agentIndices, idx1)
        else:
            agentIndices = np.append(agentIndices, idx2)
    return agentIndices

def chooseRandomMembers(populationSize):
    """
    Randomly selects two distinct indices which refer to the position of two 
    population members in the differential evolution algorithm.
    """
    idxSelection = False
    while idxSelection == False:
        idx1 = np.random.randint(0, populationSize)
        idx2 = np.random.randint(0, populationSize)
        if idx1 != idx2:
            idxSelection = True
    return idx1, idx2

def checkDETerminationCondition(numGenerations, maxGenerations, thresholdType, threshold, bestFitness, populationFitness, realComponent, imaginaryComponent, bifurType):
    """
    Checks if differential evolution has converged upon a satisfactory solution.
    If the maximum number of generations is reached or if fitness remains 
    unchanged over many iterations, the algorithm will terminate automatically. 
    The user can also specify stopping criteria for the smallest eigenvalue or 
    for the fitness value of the best member in the population by setting the
    thresholdType and threshold in evolveBifurcation. If these inputs are not specified
    by the user, evolveBifurcation will by default elect to monitor the fitness
    value to determine if this stopping criteria is satisfied.
    """
    if numGenerations == maxGenerations:
        return True
    elif (np.average(populationFitness)*1E-6 / np.std(populationFitness)) >= 1.0:
        return True
    elif thresholdType == 'eigenvalue' and checkBestEigenvals(realComponent, imaginaryComponent, bifurType) <= threshold:
        return True
    elif thresholdType == 'fitness' and bestFitness <= threshold:
        return True
    else:
        return False  
        
def checkBestEigenvals(realComponent, imaginaryComponent, bifurType):
    """
    Returns the minimum magnitude of all eigenvalue real components.
    """
    eigenvalEval = []
    if bifurType == 'oscillator':
        for i in range(len(realComponent)):
            if imaginaryComponent[i] == 0.0:
                eigenvalEval.append(10E6)
            else:
                eigenvalEval.append(abs(realComponent[i]))
    elif bifurType == 'turningPoint':
        for i in range(len(realComponent)):
            if imaginaryComponent[i] > 0.0:
                eigenvalEval.append(10E6)
            else:
                eigenvalEval.append(abs(realComponent[i]))
    if len(eigenvalEval) > 0:
        return eigenvalEval[np.argmin(eigenvalEval)]

def BFGS(objFunc, initialGuess, parameterRanges, maxIterations = 100, tol = 1E-25, BFGSArguments = (), displayBFGSProgress = True):
    """
    Performs local minimization following differential evolution using a 
    bounded Broyden-Fletcher-Goldfarb-Shanno algorithm.
    """    
    currentGuess = initialGuess
    currentFitness = objFunc(currentGuess, *BFGSArguments)
    currentHessian = np.identity(len(initialGuess))
    iterationNum = 1
    BFGSevalComplete = False
    while not BFGSevalComplete:
        currentGradient = np.dot(-1.0, approximateGradient(currentGuess, objFunc, approxGradientArguments = BFGSArguments))
        pk = np.dot(-1.0, np.dot(currentHessian, currentGradient))
        alphak = 1.0
        sk = np.dot(alphak, pk)
        updatedGuess = currentGuess + sk
        while objFunc(updatedGuess, *BFGSArguments) >= objFunc(currentGuess, *BFGSArguments):
            alphak = alphak/2.0
            sk = np.dot(alphak,pk)
            updatedGuess = currentGuess + sk
            if alphak == 0.0: 
                return currentGuess, currentFitness
        for i in range(len(updatedGuess)):
            if updatedGuess[i] < parameterRanges[i][0] or updatedGuess[i] > parameterRanges[i][1]:
                return currentGuess, currentFitness
        yk = approximateGradient(updatedGuess, objFunc, approxGradientArguments = BFGSArguments) - approximateGradient(currentGuess, objFunc, approxGradientArguments = BFGSArguments)
        updatedHessian = currentHessian + (np.dot(yk, np.transpose(yk))/np.dot(np.transpose(yk), sk)) - (np.dot(np.dot(currentHessian,sk), np.dot(np.transpose(sk), currentHessian))/np.dot(np.dot(np.transpose(sk),currentHessian), sk))
        currentGuess = updatedGuess
        currentFitness = objFunc(currentGuess, *BFGSArguments)
        currentHessian = updatedHessian
        if displayBFGSProgress == True:
            print('BFGS iteration ', str(iterationNum), ' fitness: f(x) = ', str(currentFitness))
        if currentFitness < tol or iterationNum >= maxIterations:
            BFGSevalComplete = True
        iterationNum += 1
    return currentGuess, currentFitness

def approximateGradient(parameterVector, objFunc, stepSize = 1E-8, approxGradientArguments = ()):
    """
    Computes the approximate gradient for the current point in parameter space
    by populating a vector of length len(parameterVector) with the change in
    fitness due to an incremental step divided by the size of incremental step.
    """
    startingFunctionEval = objFunc(parameterVector, *approxGradientArguments)
    step = np.ones(len(parameterVector), float)*stepSize
    gradient = np.ones(len(parameterVector), float)*((objFunc(parameterVector + step, *approxGradientArguments) - startingFunctionEval) / stepSize)
    return gradient

