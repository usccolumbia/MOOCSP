from pymoo.model.duplicate import NoDuplicateElimination
from pymoo.model.population import Population
from pymoo.model.repair import NoRepair
import numpy as np
from pymoo.model.population import pop_from_array_or_individual

class InfillCriterion:

    def __init__(self,
                 repair=None,
                 eliminate_duplicates=None,
                 n_max_iterations=100,
                 **kwargs):

        super().__init__()
        self.n_max_iterations = n_max_iterations
        self.eliminate_duplicates = eliminate_duplicates if eliminate_duplicates is not None else NoDuplicateElimination()
        self.repair = repair if repair is not None else NoRepair()
        self.parents = []
        self.f_parents=[]

    def do(self, problem, pop, n_offsprings, **kwargs):


        rand=np.random.random_sample((int(len(pop)/20),len(pop[0].get('X'))))

        off=pop_from_array_or_individual(rand)


        # infill counter - counts how often the mating needs to be done to fill up n_offsprings
        n_infills = 0
        temp=[]
        # iterate until enough offsprings are created
        while len(off) < n_offsprings:

            # how many offsprings are remaining to be created
            n_remaining = n_offsprings - len(off)

            # do the mating
            _off = self._do(problem, pop, n_remaining, **kwargs)

            _off_bk1=_off.get('X').tolist()

            # repair the individuals if necessary - disabled if repair is NoRepair
            _off = self.repair.do(problem, _off, **kwargs)

            # eliminate the duplicates - disabled if it is NoRepair
            _off = self.eliminate_duplicates.do(_off, pop, off)

            # if more offsprings than necessary - truncate them randomly
            if len(off) + len(_off) > n_offsprings:
                # IMPORTANT: Interestingly, this makes a difference in performance
                n_remaining = n_offsprings - len(off)
                _off = _off[:n_remaining]

            parents_tem=[]
            for z in _off.get('X').tolist():

                parents_tem.append(self.parents[_off_bk1.index(z)])

            temp.append(parents_tem)
            # add to the offsprings and increase the mating counter
            off = Population.merge(off, _off)
            n_infills += 1

            # if no new offsprings can be generated within a pre-specified number of generations
            if n_infills > self.n_max_iterations:
                break

        if len(temp)>0:
            temp1=[]
            for t in temp:
                temp1+=t
            self.f_parents=temp1

        return off

    def _do(self, problem, pop, n_offsprings, **kwargs):
        pass
