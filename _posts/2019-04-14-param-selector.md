---
layout: post
title: Interactive Parameter Value Comparison 
description: Comparison of our values with the literature
permalink: param_selector
---

---
The figure below can be used to explore how different parameter values change
the agreement/disagreement of the data with the theory presented in this work
and that in Razo-Mejia *et al.* 2018. Parameter values can either be entered
numerically (left column), by choosing those previously reported in the
literature (middl column), or by comparing the fits given three different values
for the allosteric energy difference. Make sure you have clicked the radio
button above the entry mode you desire before adjusting the parameters. 

The middle row of plots show the data reported in Razo-Mejia *et al.* 2018. Each
point and error corresponds to the mean and standard error of at least 10
biological replicates. The bottom left-hand plot shows leakiness measurements
presented in Garcia and Phillips 2011 (circles) and Brewster *et al.* 2014
(diamonds). The bottom right-hand plot shows gene expression measurements taken
when there are more than one regulatory architecture in the cell. These data
were collected and dissected in Brewster *et al.* 2014

This figure was generated using the [Bokeh interactive plotting
framework](http://bokeh.pydata.org). The code used to generate this figure can
be found on the [code page]({{site.baseurl}}/code/) of this website.    

<center>

{% include param_selector.html %}

</center>
