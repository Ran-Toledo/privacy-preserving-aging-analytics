# Privacy-Preserving Aging Analytics


<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
      </ul>
    </li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

![project-image]

In recent studies, it has been found that DNA methylation could be an accurate biomarker for aging in humans. The process of identifying accurate biomarkers such as this and their affecting factors could substantially increase our understanding of the aging process. Despite analytical tools to model the change of methylation and aging being relatively undeveloped, a recent line of work known as the Epigenetic PaceMaker (EPM) has developed a statistical, likelikhood-based approached to model precisely that. However, employing such a model on human subjects raises privacy concerns and if privacy is not ensured it could have a vast range of serious implications such as spikes in health insurance, job stability risk, etc. Also, the input for the EPM model relies on sensitive personal data and the model itself could also be confidential, for example; when it is a proprietary asset of a commercial entity. Therefore it is highly crucial to add a privacy component to the EPM system. In our work, we have developed privacy-preserving protocols to securely infer aging using DNAM arrays. Our solution builds on the EPM method along with privacy-preserving solutions for linear regression.


<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

This project requires the installation of the following libraries, as well as numpy, matplotlib, and scipy.

* npm
  ```sh
  npm install npm@latest -g
  ```
* phe
  ```sh
  pip install phe
  ```
* EpigeneticPacemaker
  ```sh
  pip install EpigeneticPacemaker
  ```



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

Resources used in creating this project include the following:

* [EPM Tutorial] (https://epigeneticpacemaker.readthedocs.io/en/latest/epm_tutorial/)
<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- MARKDOWN LINKS & IMAGES -->
[project-image]: https://i.im.ge/2022/08/26/OmwlhF.privacy-preserving-aging-analytics.png