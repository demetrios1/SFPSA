## Introduction and usage

### Monotone Bart

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Integration
This is the code for the integration.  Right now, you can use bart, monotone bart, or random forest, although more models will likely be added.  We also have a balanced cross-validation function for the monotone bart model, again with an option to cross validate for normal bart or random forest easily added.  

This integration takes a mean of the posterior draws and uses that as the passed probability for each obvservation in the integration optimization step.  One could also optimize the entire dataset for all the posterior draws of the joint probabilities of treatment/outcome for each observation.  However, doing so is quite the computational strain, so we do not currently include this feature. 
Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/demetrios1/Causallysensitive/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
