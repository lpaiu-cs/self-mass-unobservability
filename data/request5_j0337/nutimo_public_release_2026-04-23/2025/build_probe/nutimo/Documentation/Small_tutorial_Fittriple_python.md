<!--
SPDX-FileCopyrightText: 2024 Guillaume VOISIN (researcher at LUTH, Observatoire de Paris, PSL, CNRS) <guillaume.voisin@obspm.fr> <astro.guillaume.voisin@gmail.com>

SPDX-License-Identifier: GPL-3.0-or-later
-->
****
Here are presented some of the main functionalities of the python interface for the `Fittriple` object. It is assumed that the user is unfamiliar with Python
****

## Very basic start

* Open a python console : "ipython" is one such console, very practical and very used. 

* Import the `Fittriple` python module. In python, modules are files that contain variables, functions, object and classes. They are comparable to Fortran modules. However, the python syntax to access the content of modules is the same as that to access variables and methods of objects i.e. `module_name.function_in_module()` or `object_name.variable_in_object`.
```
    import Fittriple as ft # import the module in the file Fittriple.py under the name "ft"
```
* Create a `Fittriple` object. The name of the class in the `Fittriple` module is also `Fittriple`.
```
    fit = ft.Fittriple('TheFakePulsar/parfile-TheFakePulsar', 'TheFakePulsar.tim')  # create a "Fittriple" object from a given parfile and timfile. 
```



## Basic methods 

* In ipython, to complete a given expression one can just strike the "tabulation" key just as in a regular shell. This is very useful to know the content of the Fittriple object just created. Just type "fit." and strike "tab" to get the list of available methods.
* To get help about a particular method or object, one can write its name followed by a question mark, for instance "object.method?",  and strike "Enter". 

```
    fit.Plot_time_residuals?    # show help about "Plot_time_residuals"
    fig = fit.Plot_time_residuals() # Create a simple plot which is a "figure" object put in "fig"
    fig.show()                      # This displays the plot (not necessary to save it at the line below) 
    fig.savefig('plot.pdf')         # Saves the plot as a pdf file. Note that the method savefig recognizes the file extension given to the filename. 

    fit.Plot_time_residuals_vs_orbitphases?         # do the same with plots versus orbital phases. 
    fig = fit.Plot_time_residuals_vs_orbitphases()  
    fig.show()
    
    fit.Compute_WRMS() # Compute the weighted root-mean-square in two different ways (tempo2 computes the wrms in a way that I don't really understand. I didn't look at it too much, maybe my way is wrong).
    
    fit.Get_masses()    # Return the 3 masses (in Solar masses)
    
    timeres = fit.Get_time_residuals()  # put the time residuals in an array called "timeres"
    print timeres
    
    fit.Compute_lnposterior()   # Recompute and return the "reduced" log(Posterior probability) = -chi^2 (divided by the number of toa). Note that calling this routine is the only way to recompute all the internal quantities, in particular when a parameter is changed or a delay is enabled/disabled. 
```

## Some tools to check the validity of the computation

```
    fit.Compute_integrals_of_motion()           # Compute internally the energies and center of mass momenta (called "impulsions") and positions . Necessary before using the following lines.
    energies = fit.Get_energies()               # Get an array containing the energies
    cof_momenta = fit.Get_center_of_mass_impulsions()    # idem with momenta
    cof_positions = fit.Get_center_of_mass_positions()      #...
    fig = fit.Plot_all_conserved_quantities()               # Get a plot 
    fig.show()                                              # and show it (it could be saved as well). 
    
    bats, delay_geom, delay_ein, delay_shap, delay_aber = fit.Get_fake_bats_and_delays_interp(10000) # Generate 10000 fake bats evenly spaced on the current time span assuming the current parameters and return as well the corresponding delays (depending on the enabled delays). 
    fig_geometric, fig_einstein, fig_shapiro, fig_aberration = fit.Plot_delays() # Create four figures of the residuals minus a given time delay. 
```

** NB : the previous line could be written "figs = fit.Plot_delays()". In this case "figs" is a python tuple that contains four objects (that could of heterogeneous type). The elements of the tuple can then be accessed with an array-like syntax : for example "fig_einstein = figs[1]" (array indexes start at 0 in python). This applies to any python function with multiple output. 
**

## Enabling / disabling effects 
```
    fit.Get_theory() # get the theory used for the equations of motion : if zero is returned then the equations are Newtonian, if 1 is returned then 1PN equations are used.
    fit.Set_theory(1) # Set the theory to post-Newtonian
    
    fit.Set_shapiro_delay_off() # disable the Shapiro delay
    fit.Compute_lnposterior()   # Recompute with Shapiro disabled
    fig = fit.Plot_time_residuals()     # Plot the time delays 
    fig.show()  
    fit.Set_shapiro_delay_on()  # Re-enable the Shapiro delay
    fit.Compute_lnposterior()   # Recompute
```

**NB : The other delays can be enabled/disabled with similar methods.**