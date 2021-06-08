class Node:
	def __init__(self,aaid,pos):
		self.id = [aaid]
		self.position = [pos]
	def __hash__(self):
		return hash((str(self.id),str(self.position)))
	def __eq__(self, other):
		if not isinstance(other, type(self)):
			return NotImplemented
		return self.id == other.id and self.position == other.position
	def __lt__(self, other):
		if not isinstance(other, type(self)):
			return NotImplemented
		return self.position < other.position
	def __str__(self):
		return properties[self.id[0]] + str(self.position[0])
	def add2Cluster(self,residue):
		self.id.append(residue)
	def getClusters(self):
		cluster = []
		for i in range(0,len(self.id)):
			cluster.append(properties[self.id[i]] + str(self.position[i]))
		return cluster
	def toString(self):
		if len(self.id) > 1:
			return str(self.id)
		return properties[self.id[0]] + str(self.position[0])
	def toStackedString(self):
		if len(self.id) > 1:
			return str(self.id)
		return properties[self.id[0]] + "_" + str(self.position[0])
	def getTuple(self):
		return (self.id[0],self.position[0])
	def getIndex(self):
		return self.position[0]-1