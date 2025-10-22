"""
Definition of _attribute and _attribute_collection classes.
"""

from unduwave.unduwave_incl import *

class _attribute:
	"""
	Represents an attribute with a value.

	Args:
		value: Initial value of the attribute.
		name (str, optional): Name of the attribute.
	"""
	def __init__(self, value=None, name=None, in_name = None,unit='',fac=None):
		"""
		Initialize an attribute.
		:param value: The value held by this class
		:param str name: Name of the attribute
		:param str in_name: Name as used by Wave or Undumag
		:param str unit: The physical unit of the quantity
		'param float fac': gauging factor
		"""
		self._value = value
		self.set_name(name)
		self._in_name = in_name
		self._unit=unit
		self._fac=fac

	def set(self, value):
		"""
		Sets the value
		"""
		self._value = value

	def get(self):
		"""
		Returns the value
		"""
		return self._value

	def get_fac(self):
		"""
		Returns scaling factor
		"""
		if self._fac is None :
			return 1
		return self._fac

	def set_name(self,name) :
		"""
		Setting the name
		"""
		self._name = name

	def set_in_name(self, in_name):
		"""
		Setting the in-name
		"""
		self._in_name = in_name

	def get_in_name(self):
		"""
		Returns in-name
		"""
		return self._in_name

	@property
	def name(self):
		return self._name

	def __str__(self):
		return str(self._value)

	def __repr__(self):
		return f"{self.__class__.__name__}('{self}')"

	def __iter__(self):
		for attr, value in self.__dict__.items():
			yield attr, value

class _attribute_collection:
	"""
	Represents a collection of attributes with children.
	This is a class containing attributes as class-members.
	The class members are administered in unison.
	"""
	def __init__(self):
		"""
		Initialize the attribute collection
		"""
		self._children = []
		self._add_attributes()

	def _add_attributes(self):
		"""
		Scans all class members, finds those of type _attribute and sets those 
		as new attributes
		"""

		for attr_name, attr_value in dict(self).items():
			if isinstance(attr_value, _attribute):
				attr_value.set_name(attr_name) # set the name of the copy
				self._children.append(getattr(self, attr_name))

		# attributes = vars(self.__class__)
		# for attr_name, attr_value in attributes.items():
		# 	if isinstance(attr_value, _attribute):
		# 		attr = copy.copy(attr_value) # copy the attribute with all its attributes
		# 		attr.set_name(attr_name) # set the name of the copy
		# 		setattr(self, attr_name, attr) # set the copied version as the class attribute - auto setting of name attribute
		# 		# setattr(self, attr_name, WaveAttribute(attr_value.get(), name=attr_name, wave_in_name=attr_value.get_wave_in_name()))
		# 		self._children.append(getattr(self, attr_name))

	def show_all_children(self) : 
		"""
		Prints the names of the children
		"""
		for child in self.children() :
			print(f"Para {child.get_in_name()}={child.get()}")

	def children(self):
		"""
		Yields the children of the _attribute_collection.
		"""
		for child in self._children:
			yield child

	def __iter__(self):
		for attr, value in self.__dict__.items():
			yield attr, value
