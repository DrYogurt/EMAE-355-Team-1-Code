# The Fluid Class

I've created a wrapper for the CoolProp class in order to more easily create and find the state of fluids. This wrapper has two benefits:

1. it uses pint for units, which limits errors in arithmetic
2. it has a nicer notation, using the `__get_attr__` magic keyword.

For any fluid, simply call `fluid.property(state1=val1,state2=val2)`, this will return the correct value for that property, for instance:

```python
co2 = Fluid('co2')
d = co2.density(T=300, P=1e5)
print(d)
```

prints:

```
1.773026407279758 kilogram / meter ** 3
```