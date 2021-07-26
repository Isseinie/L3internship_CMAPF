class PriorityQueue:

    def __init__(self, nb_elements_max):
        self.heapsize = 0
        self.heap = [(None, None) for i in range(nb_elements_max)]
        self.dict = [-1 for element in range(nb_elements_max)]
    
    def pop(self):
        min = self.heap[0]
        self.dict[self.heap[self.heapsize-1][0]]=0
        self.heap[0] = self.heap[self.heapsize-1]
        self.heapsize-=1
        self._bubbledown(0)
        return min[0]

    def _bubbledown(self, i):
        j = i
        if 2*i+1< self.heapsize:
            if self.heap[2*i+1][1]<self.heap[j][1]:
                j = 2*i+1
        if 2*i+2< self.heapsize:
            if self.heap[2*i+2][1]<self.heap[j][1]:
                j = 2*i+2
        if i!=j:
            self._swap(i, j)
            self._bubbledown(j)
    
    def _swap(self, i, j):
        self.heap[i], self.heap[j] = self.heap[j], self.heap[i]
        self.dict[self.heap[i][0]] = i
        self.dict[self.heap[j][0]] = j

    def push(self, element, priority):
        self.heapsize+=1
        self.heap[self.heapsize-1]=[element, priority]
        self.dict[element] = self.heapsize-1
        self._bubbleup(self.heapsize-1)
    
    def _bubbleup(self, i):
        if i >0:
            parent = (i-1)//2
            if self.heap[i][1] < self.heap[parent][1]:
                self._swap(i, parent)
                self._bubbleup(parent)
    
    def decrease(self, element, priority):
        if self.dict[element]==-1:
            self.push(element, priority)
        else :
            self.heap[self.dict[element]][1]= priority
            self._bubbleup(self.dict[element])

    def is_empty(self):
        return (self.heapsize==0)