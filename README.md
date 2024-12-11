# From Big Data to Small Scales: Machine Learning Enhances Microclimate Model Predictions

_Alon Itzkovitch, Idan Sulami, Ronny Doron Efroni, Moni Shahar, and
Ofir Levy_

### Please contact Ofir Levy (levyofir@tauex.tau.ac.il) about the code or data

## Abstract:

1) Microclimates are critical for understanding how organisms interact with their environments, shaping behaviour, physiology, and species distributions. They influence key ecological processes, including thermoregulation, habitat use, and population dynamics. Traditional approaches to modelling ground temperatures in microhabitats primarily rely on physical heat-balance models informed by remote sensing data. While these models have been instrumental in ecological research, they often exhibit biases due to unaccounted environmental complexities, poorly constrained parameters, and inherent simplifications in the physical processes.

2) In this study, we introduce a novel framework that integrates high-resolution drone-based microclimate mapping with machine learning to improve microclimate model predictions. Using drone imagery at a spatial resolution of 15 cm, we generated detailed environmental maps, including solar radiation, vegetation indices, and skyview factors, to parameterize a physical heat-balance microclimate model for predicting ground temperatures. We validated the model’s predictions using thermal maps derived from a drone-mounted infrared camera. To address model errors, we trained a machine learning model, based on a random-forest algorithm, to predict and reduce errors for new prediction maps. 

3) Our machine learning model effectively reduced systematic errors (biases) and improved prediction accuracy across diverse environmental conditions. We identified key factors driving microclimate model inaccuracies, such as vegetation cover, solar radiation, and height above ground, offering actionable insights for refining future physical models. The machine learning correction reduced mean errors by over 30\% and consistently narrowed the range of prediction errors, demonstrating robust performance across test datasets.

4) By integrating machine learning with traditional physical modelling, we bridge a long-standing gap in microclimate research. Our approach advances the field into the era of big data, providing a scalable tool for ecologists, conservation practitioners, and land managers. This framework enables the generation of accurate, fine-scale microclimate predictions that are essential for understanding species responses to climate change, designing habitat conservation strategies, and informing climate-resilient management practices in a rapidly changing world. \end{enumerate}


# **Repository Directory**:
## See the `Example data` subdirectory for data and metadata. <!-- [link](https://github.com/levyofi/Itzkovitch_el_al_Proceedings_B/tree/main/Example%20data). -->

## See the `Code` subdirectory for codes.  <!--: [link](https://github.com/levyofi/Itzkovitch_el_al_Proceedings_B/tree/main/Codes).-->
