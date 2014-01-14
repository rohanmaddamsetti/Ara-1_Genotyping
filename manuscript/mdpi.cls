%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MDPI      class for LaTeX files              21-04-2011
%% %%
%% %% class style provided by Stefano Mariani
%% %% for any information please send an e-mail to:
%% %%                  stefano.mariani@polimi.it
%% %%   Modified by Dietrich Rofdorf (rordorf@mdpi.com)
%% %%              Peter Harremoes (entropeter@mdpi.com)
%% %%              Zeno Schumacher (latex@mdpi.com)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{mdpi}[2011/01/13 MDPI paper class]

\AtEndDocument{\par \cright}
% % Choose journal
\DeclareOption{administrativesciences}{\gdef\@journal{AdministrativeSciences}}
\DeclareOption{algorithms}{\gdef\@journal{Algorithms}}
\DeclareOption{atmosphere}{\gdef\@journal{Atmosphere}}
\DeclareOption{cancers}{\gdef\@journal{Cancers}}
\DeclareOption{challenges}{\gdef\@journal{Challenges}}
\DeclareOption{computers}{\gdef\@journal{Computers}}
\DeclareOption{diagnostics}{\gdef\@journal{Diagnostics}}
\DeclareOption{diversity}{\gdef\@journal{Diversity}}
\DeclareOption{energies}{\gdef\@journal{Energies}}
\DeclareOption{entropy}{\gdef\@journal{Entropy}}
\DeclareOption{forests}{\gdef\@journal{Forests}}
\DeclareOption{futureinternet}{\gdef\@journal{FutureInternet}}
\DeclareOption{games}{\gdef\@journal{Games}}
\DeclareOption{genes}{\gdef\@journal{Genes}}
\DeclareOption{information}{\gdef\@journal{Information}}
\DeclareOption{ijerph}{\gdef\@journal{Ijerph}}
\DeclareOption{ijms}{\gdef\@journal{Ijms}}
\DeclareOption{jlpea}{\gdef\@journal{Jlpea}}
\DeclareOption{marinedrugs}{\gdef\@journal{MarineDrugs}}
\DeclareOption{materials}{\gdef\@journal{Materials}}
\DeclareOption{micromachines}{\gdef\@journal{Micromachines}}
\DeclareOption{molbank}{\gdef\@journal{Molbank}}
\DeclareOption{molecules}{\gdef\@journal{Molecules}}
\DeclareOption{nutrients}{\gdef\@journal{Nutrients}}
\DeclareOption{pharmaceuticals}{\gdef\@journal{Pharmaceuticals}}
\DeclareOption{pharmaceutics}{\gdef\@journal{Pharmaceutics}}
\DeclareOption{polymers}{\gdef\@journal{Polymers}}
\DeclareOption{religions}{\gdef\@journal{Religions}}
\DeclareOption{remotesensing}{\gdef\@journal{RemoteSensing}}
\DeclareOption{sensors}{\gdef\@journal{Sensors}}
\DeclareOption{sustainability}{\gdef\@journal{Sustainability}}
\DeclareOption{symmetry}{\gdef\@journal{Symmetry}}
\DeclareOption{toxins}{\gdef\@journal{Toxins}}
\DeclareOption{viruses}{\gdef\@journal{Viruses}}
\DeclareOption{water}{\gdef\@journal{Water}}
%%%%Uppercase
%%%%will have the same effect
\DeclareOption{Administrativesciences}{\gdef\@journal{AdministrativeSciences}}
\DeclareOption{Algorithms}{\gdef\@journal{Algorithms}}
\DeclareOption{Atmosphere}{\gdef\@journal{Atmosphere}}
\DeclareOption{Cancers}{\gdef\@journal{Cancers}}
\DeclareOption{Challenges}{\gdef\@journal{Challenges}}
\DeclareOption{Computers}{\gdef\@journal{Computers}}
\DeclareOption{Diagnostics}{\gdef\@journal{Diagnostics}}
\DeclareOption{Diversity}{\gdef\@journal{Diversity}}
\DeclareOption{Energies}{\gdef\@journal{Energies}}
\DeclareOption{Entropy}{\gdef\@journal{Entropy}}
\DeclareOption{Forests}{\gdef\@journal{Forests}}
\DeclareOption{Futureinternet}{\gdef\@journal{FutureInternet}}
\DeclareOption{Games}{\gdef\@journal{Games}}
\DeclareOption{Genes}{\gdef\@journal{Genes}}
\DeclareOption{Information}{\gdef\@journal{Information}}
\DeclareOption{Ijerph}{\gdef\@journal{Ijerph}}
\DeclareOption{Ijms}{\gdef\@journal{Ijms}}
\DeclareOption{Jlpea}{\gdef\@journal{Jlpea}}
\DeclareOption{MarineDrugs}{\gdef\@journal{MarineDrugs}}
\DeclareOption{Materials}{\gdef\@journal{Materials}}
\DeclareOption{Micromachines}{\gdef\@journal{Micromachines}}
\DeclareOption{Molbank}{\gdef\@journal{Molbank}}
\DeclareOption{Molecules}{\gdef\@journal{Molecules}}
\DeclareOption{Nutrients}{\gdef\@journal{Nutrients}}
\DeclareOption{Pharmaceuticals}{\gdef\@journal{Pharmaceuticals}}
\DeclareOption{Pharmaceutics}{\gdef\@journal{Pharmaceutics}}
\DeclareOption{Polymers}{\gdef\@journal{Polymers}}
\DeclareOption{Religions}{\gdef\@journal{Religions}}
\DeclareOption{RemoteSensing}{\gdef\@journal{RemoteSensing}}
\DeclareOption{Sensors}{\gdef\@journal{Sensors}}
\DeclareOption{Sustainability}{\gdef\@journal{Sustainability}}
\DeclareOption{Symmetry}{\gdef\@journal{Symmetry}}
\DeclareOption{Toxins}{\gdef\@journal{Toxins}}
\DeclareOption{Viruses}{\gdef\@journal{Viruses}}
\DeclareOption{Water}{\gdef\@journal{Water}}
\DeclareOption{journal}{\ClassWarning{mdpi}{You have not
specified the journal. The first class option of the document
class mdpi specify the journal. The word journal should be
replaced by one of the following journal names: administrativesciences, algorithms, atmosphere, cancers, challenges,
diversity, energies, entropy, forests, futureinternet, games, genes, information, ijerph, ijms, marinedrugs,
materials, micromachines, molbank, molecules, nutrients, pharmaceuticals, pharmaceutics , polymers, religions, remotesensing,
sensors, sustainability, symmetry, toxins, viruses, water}}
% % Choose type of manuscript
\DeclareOption{article}{\gdef\@arttype{Article}}
\DeclareOption{review}{\gdef\@arttype{Review}}
\DeclareOption{communication}{\gdef\@arttype{Communication}}
\DeclareOption{commentary}{\gdef\@arttype{Commentary}}
\DeclareOption{bookreview}{\gdef\@arttype{Book Review}}
\DeclareOption{correction}{\gdef\@arttype{Correction}}
\DeclareOption{addendum}{\gdef\@arttype{Addendum}}
\DeclareOption{editorial}{\gdef\@arttype{Editorial}}
\DeclareOption{changes}{\gdef\@arttype{Changes}}

% Choose the status of the manuscript
\DeclareOption{submit}{\gdef\@status{submit}}
\DeclareOption{accept}{\gdef\@status{accept}}

% Choose the whether there is one or more authors
\DeclareOption{oneauthor}{\gdef\@authornum{author}}
\DeclareOption{moreauthors}{\gdef\@authornum{authors}}
\DeclareOption{10pt}{\OptionNotUsed}
\DeclareOption{twocolumn}{\OptionNotUsed}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}

\ExecuteOptions{entropy,12pt,a4paper,article,submit,oneauthor}
\ProcessOptions\relax

\LoadClass[12pt]{article}

\RequirePackage{calc}
\RequirePackage{indentfirst}
\RequirePackage{times,fancyhdr}
\RequirePackage{graphicx}
\RequirePackage{soul}
\RequirePackage{lastpage}
\RequirePackage{ifthen}
\RequirePackage{lineno}
\RequirePackage{float}
\floatstyle{plaintop}
\restylefloat{figure}
\restylefloat{table}
\usepackage{caption}
\captionsetup{position=top, margin=1cm, labelfont=bf, labelsep=period, skip=12pt, justification=justified}
\RequirePackage{color}
\ifthenelse{\equal{\@journal}{AdministrativeSciences}}{\definecolor{logocolor}{rgb}{0.31,0.51,0.60}}{}%
\ifthenelse{\equal{\@journal}{Algorithms}}{\definecolor{logocolor}{rgb}{0.192,0.383,0}}{}%
\ifthenelse{\equal{\@journal}{Atmosphere}}{\definecolor{logocolor}{rgb}{0.2,0.2,0.6}}{}%
\ifthenelse{\equal{\@journal}{Cancers}}{\definecolor{logocolor}{rgb}{0,0.20,0.60}}{}%
\ifthenelse{\equal{\@journal}{Challenges}}{\definecolor{logocolor}{rgb}{0.13,0.44,0.23}}{}%
\ifthenelse{\equal{\@journal}{Computers}}{\definecolor{logocolor}{rgb}{0,0,0}}{}%
\ifthenelse{\equal{\@journal}{Diagnostics}}{\definecolor{logocolor}{rgb}{0.17,0.51,0.56}}{}%
\ifthenelse{\equal{\@journal}{Diversity}}{\definecolor{logocolor}{rgb}{0.25,0.41,0.49}}{}%
\ifthenelse{\equal{\@journal}{Energies}}{\definecolor{logocolor}{rgb}{0.6,0.2,0}}{}%
\ifthenelse{\equal{\@journal}{Entropy}}{\definecolor{logocolor}{rgb}{0,0.502,0.502}}{}%
\ifthenelse{\equal{\@journal}{Forests}}{\definecolor{logocolor}{rgb}{0.16,0.31,0.13}}{}%
\ifthenelse{\equal{\@journal}{FutureInternet}}{\definecolor{logocolor}{rgb}{0.35,0.35,0.35}}{}%
\ifthenelse{\equal{\@journal}{Games}}{\definecolor{logocolor}{rgb}{0.34,0.18,0.14}}{}%
\ifthenelse{\equal{\@journal}{Genes}}{\definecolor{logocolor}{rgb}{0.37,0.29,0.48}}{}%
\ifthenelse{\equal{\@journal}{Information}}{\definecolor{logocolor}{rgb}{0.03,0.5,0.27}}{}%
\ifthenelse{\equal{\@journal}{Ijerph}}{\definecolor{logocolor}{rgb}{0,0.502,0}}{}%
\ifthenelse{\equal{\@journal}{Ijms}}{\definecolor{logocolor}{rgb}{0,0.502,0.502}}{}%
\ifthenelse{\equal{\@journal}{Jlpea}}{\definecolor{logocolor}{rgb}{0.07,0.37,0.22}}{}%
\ifthenelse{\equal{\@journal}{MarineDrugs}}{\definecolor{logocolor}{rgb}{0,0.502,0.502}}{}%
\ifthenelse{\equal{\@journal}{Materials}}{\definecolor{logocolor}{rgb}{0.6,0.4,0.2}}{}%
\ifthenelse{\equal{\@journal}{Micromachines}}{\definecolor{logocolor}{rgb}{0,0.2,0.6}}{}%
\ifthenelse{\equal{\@journal}{Molecules}}{\definecolor{logocolor}{rgb}{0,0.502,0.502}}{}%
\ifthenelse{\equal{\@journal}{Molbank}}{\definecolor{logocolor}{rgb}{0,0.502,0.502}}{}%
\ifthenelse{\equal{\@journal}{Nutrients}}{\definecolor{logocolor}{rgb}{0.20,0.40,0.20}}{}%
\ifthenelse{\equal{\@journal}{Pharmaceuticals}}{\definecolor{logocolor}{rgb}{0.44,0.19,0.63}}{}%
\ifthenelse{\equal{\@journal}{Pharmaceutics}}{\definecolor{logocolor}{rgb}{0.20,0.40,0.40}}{}%
\ifthenelse{\equal{\@journal}{Polymers}}{\definecolor{logocolor}{rgb}{0.44,0.19,0.63}}{}%
\ifthenelse{\equal{\@journal}{Religions}}{\definecolor{logocolor}{rgb}{0,0.4,0.6}}{}%
\ifthenelse{\equal{\@journal}{RemoteSensing}}{\definecolor{logocolor}{rgb}{0.20,0.40,0}}{}%
\ifthenelse{\equal{\@journal}{Sensors}}{\definecolor{logocolor}{rgb}{0,0.502,0.502}}{}%
\ifthenelse{\equal{\@journal}{Sustainability}}{\definecolor{logocolor}{rgb}{0.21,0.37,0.57}}{}%
\ifthenelse{\equal{\@journal}{Symmetry}}{\definecolor{logocolor}{rgb}{0,0.20,0.20}}{}%
\ifthenelse{\equal{\@journal}{Toxins}}{\definecolor{logocolor}{rgb}{0.60,0.20,0}}{}%
\ifthenelse{\equal{\@journal}{Viruses}}{\definecolor{logocolor}{rgb}{0.5,0,0}}{}%
\ifthenelse{\equal{\@journal}{Water}}{\definecolor{logocolor}{rgb}{0.21,0.37,0.57}}{}%
%\definecolor{greens}{rgb}{0,0.502,0.502}
\RequirePackage[sort&compress]{natbib}
\bibpunct{[}{]}{,}{n}{}{}

\renewcommand\NAT@set@cites{%
  \ifNAT@numbers
    \ifNAT@super \let\@cite\NAT@citesuper
       \def\NAT@mbox##1{\unskip\nobreak\textsuperscript{##1}}%
       \let\citeyearpar=\citeyear
       \let\NAT@space\relax
       \def\NAT@super@kern{\kern\p@}%
    \else
       \let\NAT@mbox=\mbox
       \let\@cite\NAT@citenum
       \let\NAT@space\relax
       \let\NAT@super@kern\relax
    \fi
    \let\@citex\NAT@citexnum
    \let\@biblabel\NAT@biblabelnum
    \let\@bibsetup\NAT@bibsetnum
    \renewcommand\NAT@idxtxt{\NAT@name\NAT@spacechar\NAT@open\NAT@num\NAT@close}%
    \def\natexlab##1{}%
    \def\NAT@penalty{\penalty\@m}%
  \else
    \let\@cite\NAT@cite
    \let\@citex\NAT@citex
    \let\@biblabel\NAT@biblabel
    \let\@bibsetup\NAT@bibsetup
    \let\NAT@space\NAT@spacechar
    \let\NAT@penalty\@empty
    \renewcommand\NAT@idxtxt{\NAT@name\NAT@spacechar\NAT@open\NAT@date\NAT@close}%
    \def\natexlab##1{##1}%
  \fi}

\renewcommand{\baselinestretch}{1.19}

%%%%%%%%%
%%
%%hyperref package and options
%%%%%%%%

\RequirePackage[unicode=false, pdftoolbar=true, pdfmenubar=true, bookmarksopen={true},
pdffitwindow=true, pdfnewwindow=true, colorlinks=true,
linkcolor=blue, citecolor=blue, urlcolor=blue]{hyperref}
\newcommand{\org@maketitle}{}% LATEX-Check
\let\org@maketitle\maketitle
\def\maketitle{%
\hypersetup{
pdftitle={\@Title},
pdfsubject={\@abstract},
pdfkeywords={\@keyword}
}%
\org@maketitle
}
\renewcommand{\@pdfauthor}{\@Author}

\RequirePackage{url}
\urlstyle{same}
\clubpenalty=10000
\widowpenalty=10000
\displaywidowpenalty=10000

%%%%For tracking changes
\ifthenelse{\equal{\@status}{accept}}{\newcommand{\cor}[2]{{#2}}}{
\newcommand{\cor}[2]{{\textcolor{red}{\st{#1 }}}\textcolor[rgb]{0.25,0,0.5}{#2}}}

%%%

\newcommand{\firstargument}{}
\newcommand{\Title}[1]{\gdef\@Title{#1}}%
\newcommand{\Author}[1]{\gdef\@Author{#1}}%
\newcommand{\lastpage}[1]{\gdef\@lastpage{#1}}
\newcommand{\doinum}[1]{\gdef\@doinum{#1}}
\newcommand{\pubvolume}[1]{\gdef\@pubvolume{#1}}
\newcommand{\pubyear}[1]{\gdef\@pubyear{#1}}
\newcommand{\address}[2][]{\renewcommand{\firstargument}{#1}\gdef\@address{#2}}
\newcommand{\corresfirstargument}{}
\newcommand{\corres}[2][]{\renewcommand{\corresfirstargument}{#1}\gdef\@corres{#2}}
\renewcommand{\abstract}[1]{\gdef\@abstract{#1}}
\newcommand{\history}[1]{\gdef\@history{#1}} \def\@pacs{}
\newcommand{\PACS}[1]{\gdef\@pacs{#1}} \def\@msc{}
\newcommand{\MSC}[1]{\gdef\@msc{#1}}
\newcommand{\keyword}[1]{\gdef\@keyword{#1}}

%% Define the journal name for the header using the ifthen
%%package
\newcommand{\journalname}{%
  \ifthenelse{\equal{\@journal}{AdministrativeSciences}}{Adm. Sci.}{}%
  \ifthenelse{\equal{\@journal}{Algorithms}}{Algorithms}{}%
  \ifthenelse{\equal{\@journal}{Atmosphere}}{Atmosphere}{}%
  \ifthenelse{\equal{\@journal}{Cancers}}{Cancers}{}%
  \ifthenelse{\equal{\@journal}{Challenges}}{Challenges}{}%
  \ifthenelse{\equal{\@journal}{Computers}}{Computers}{}%
  \ifthenelse{\equal{\@journal}{Diagnostics}}{Diagnostics}{}%
  \ifthenelse{\equal{\@journal}{Diversity}}{Diversity}{}%
  \ifthenelse{\equal{\@journal}{Energies}}{Energies}{}%
  \ifthenelse{\equal{\@journal}{Entropy}}{Entropy}{}%
  \ifthenelse{\equal{\@journal}{Forests}}{Forests}{}%
  \ifthenelse{\equal{\@journal}{FutureInternet}}{Future Internet}{}%
  \ifthenelse{\equal{\@journal}{Games}}{Games}{}%
  \ifthenelse{\equal{\@journal}{Genes}}{Genes}{}%
  \ifthenelse{\equal{\@journal}{Information}}{Information}{}%
  \ifthenelse{\equal{\@journal}{Ijerph}}{Int. J. Environ. Res. Public Health}{}%
  \ifthenelse{\equal{\@journal}{Ijms}}{Int. J. Mol. Sci.}{}%
  \ifthenelse{\equal{\@journal}{Jlpea}}{J. Low Power Electron. Appl.}{}%
  \ifthenelse{\equal{\@journal}{MarineDrugs}}{Mar. Drugs}{}%
  \ifthenelse{\equal{\@journal}{Materials}}{Materials}{}%
  \ifthenelse{\equal{\@journal}{Micromachines}}{Micromachines}{}%
  \ifthenelse{\equal{\@journal}{Molbank}}{Molbank}{}%
  \ifthenelse{\equal{\@journal}{Molecules}}{Molecules}{}%
  \ifthenelse{\equal{\@journal}{Nutrients}}{Nutriens}{}%
  \ifthenelse{\equal{\@journal}{Pharmaceuticals}}{Pharmaceuticals}{}%
  \ifthenelse{\equal{\@journal}{Pharmaceutics}}{Pharmaceutics}{}%
  \ifthenelse{\equal{\@journal}{Polymers}}{Polymers}{}%
  \ifthenelse{\equal{\@journal}{Religions}}{Religions}{}%
  \ifthenelse{\equal{\@journal}{RemoteSensing}}{Remote Sens.}{}%
  \ifthenelse{\equal{\@journal}{Sensors}}{Sensors}{}%
  \ifthenelse{\equal{\@journal}{Sustainability}}{Sustainability}{}%
  \ifthenelse{\equal{\@journal}{Symmetry}}{Symmetry}{}%
  \ifthenelse{\equal{\@journal}{Toxins}}{Toxins}{}%
  \ifthenelse{\equal{\@journal}{Viruses}}{Viruses}{}%
  \ifthenelse{\equal{\@journal}{Water}}{Water}{}%
  }
%% Define the logo size of the journal using the ifthen
%%%package
\newcommand{\logosize}{%
  \ifthenelse{\equal{\@journal}{AdministrativeSciences}}{\fontsize{20}{20}}{}%
  \ifthenelse{\equal{\@journal}{Algorithms}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Atmosphere}}{\fontsize{30}{30}}{}%
  \ifthenelse{\equal{\@journal}{Cancers}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Challenges}}{\fontsize{36}{36}}{}%
  \ifthenelse{\equal{\@journal}{Computers}}{\fontsize{28}{28}}{}%
  \ifthenelse{\equal{\@journal}{Diagnostics}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Diversity}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Energies}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Entropy}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Forests}}{\fontsize{36}{36}}{}%
  \ifthenelse{\equal{\@journal}{FutureInternet}}{\fontsize{28}{28}}{}%
  \ifthenelse{\equal{\@journal}{Games}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Genes}}{\fontsize{36}{36}}{}%
  \ifthenelse{\equal{\@journal}{Information}}{\fontsize{28}{28}}{}%
  \ifthenelse{\equal{\@journal}{Ijerph}}{\fontsize{14}{14}}{}%
  \ifthenelse{\equal{\@journal}{Ijms}}{\fontsize{14}{14}}{}%
  \ifthenelse{\equal{\@journal}{Jlpea}}{\fontsize{14}{14}}{}%
  \ifthenelse{\equal{\@journal}{MarineDrugs}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Materials}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Micromachines}}{\fontsize{26}{26}}{}%
  \ifthenelse{\equal{\@journal}{Molecules}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Molbank}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Nutrients}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Pharmaceuticals}}{\fontsize{26}{26}}{}%
  \ifthenelse{\equal{\@journal}{Pharmaceutics}}{\fontsize{26}{26}}{}%
  \ifthenelse{\equal{\@journal}{Polymers}}{\fontsize{36}{36}}{}%
  \ifthenelse{\equal{\@journal}{Religions}}{\fontsize{36}{36}}{}%
  \ifthenelse{\equal{\@journal}{RemoteSensing}}{\fontsize{24}{24}}{}%
  \ifthenelse{\equal{\@journal}{Sensors}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Sustainability}}{\fontsize{28}{28}}{}%
  \ifthenelse{\equal{\@journal}{Symmetry}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Toxins}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Viruses}}{\fontsize{32}{32}}{}%
  \ifthenelse{\equal{\@journal}{Water}}{\fontsize{36}{36}}{}%
  }
%% % % Creates the command that makes the logo of the journal.
\newcommand{\logo}{%
    \ifthenelse{\equal{\@journal}{AdministrativeSciences}}{\\[-32pt]\fontsize{20}{20}\bfseries\itshape{administrative}\\[4pt] \bfseries\itshape\fontsize{35}{35}\selectfont{sciences}}{}%
    \ifthenelse{\equal{\@journal}{Algorithms}}{algorithms}{}%
    \ifthenelse{\equal{\@journal}{Atmosphere}}{atmosphere}{}%
    \ifthenelse{\equal{\@journal}{Cancers}}{cancers}{}%
    \ifthenelse{\equal{\@journal}{Challenges}}{challenges}{}%
    \ifthenelse{\equal{\@journal}{Computers}}{computers}{}%
    \ifthenelse{\equal{\@journal}{Diagnostics}}{diagnostics}{}%
    \ifthenelse{\equal{\@journal}{Diversity}}{diversity}{}%
    \ifthenelse{\equal{\@journal}{Energies}}{energies}{}%
    \ifthenelse{\equal{\@journal}{Entropy}}{entropy}{}%
    \ifthenelse{\equal{\@journal}{Forests}}{forests}{}%
    \ifthenelse{\equal{\@journal}{FutureInternet}}{future internet}{}%
    \ifthenelse{\equal{\@journal}{Games}}{games}{}%
    \ifthenelse{\equal{\@journal}{Genes}}{genes}{}%
    \ifthenelse{\equal{\@journal}{Information}}{information}{}%
    \ifthenelse{\equal{\@journal}{Ijerph}}{\normalfont{International Journal of}\\\normalfont\bfseries\fontsize{16}{16}{Environmental Research and\\ Public Health}}{}%
    \ifthenelse{\equal{\@journal}{Ijms}}{\normalfont\fontsize{13}{13}{International Journal of}\\\normalfont\bfseries\fontsize{18}{18}\selectfont{Molecular Sciences}}{}%
    \ifthenelse{\equal{\@journal}{Jlpea}}{\normalfont\bfseries\fontsize{14}{14}{Journal of}\\\normalfont\bfseries\fontsize{14}{14}\selectfont{Low Power Electronics}\\\normalfont\bfseries\fontsize{14}{14}\selectfont{and Applications}}{}%
    \ifthenelse{\equal{\@journal}{MarineDrugs}}{Marine Drugs}{}%
    \ifthenelse{\equal{\@journal}{Materials}}{materials}{}%
    \ifthenelse{\equal{\@journal}{Micromachines}}{micromachines}{}%
    \ifthenelse{\equal{\@journal}{Molecules}}{molecules}{}%
    \ifthenelse{\equal{\@journal}{Molbank}}{molbank}{}%
    \ifthenelse{\equal{\@journal}{Nutrients}}{nutrients}{}%
    \ifthenelse{\equal{\@journal}{Pharmaceuticals}}{pharmaceuticals}{}%
    \ifthenelse{\equal{\@journal}{Pharmaceutics}}{pharmaceutics}{}%
    \ifthenelse{\equal{\@journal}{Polymers}}{polymers}{}%
    \ifthenelse{\equal{\@journal}{Religions}}{religions}{}%
    \ifthenelse{\equal{\@journal}{RemoteSensing}}{Remote Sensing}{}%
    \ifthenelse{\equal{\@journal}{Sensors}}{sensors}{}%
    \ifthenelse{\equal{\@journal}{Sustainability}}{sustainability}{}%
    \ifthenelse{\equal{\@journal}{Symmetry}}{symmetry}{}%
    \ifthenelse{\equal{\@journal}{Toxins}}{toxins}{}%
    \ifthenelse{\equal{\@journal}{Viruses}}{viruses}{}%
    \ifthenelse{\equal{\@journal}{Water}}{water}{}%
}
% % Creates the ISSN for the journal using the ifthen package
\newcommand{\issn}{%
  \ifthenelse{\equal{\@journal}{AdministrativeSciences}}{2076-3387}{}%
  \ifthenelse{\equal{\@journal}{Algorithms}}{1999-4893}{}%
  \ifthenelse{\equal{\@journal}{Atmosphere}}{2073-4433}{}%
  \ifthenelse{\equal{\@journal}{Cancers}}{2072-6694}{}%
  \ifthenelse{\equal{\@journal}{Challenges}}{2078-1547}{}%
  \ifthenelse{\equal{\@journal}{Computers}}{2073-431X}{}%
  \ifthenelse{\equal{\@journal}{Diagnostics}}{2075-4418}{}%
  \ifthenelse{\equal{\@journal}{Diversity}}{1424-2818}{}%
  \ifthenelse{\equal{\@journal}{Energies}}{1996-1073}{}%
  \ifthenelse{\equal{\@journal}{Entropy}}{1099-4300}{}%
  \ifthenelse{\equal{\@journal}{Forests}}{1999-4907}{}%
  \ifthenelse{\equal{\@journal}{FutureInternet}}{1999-5903}{}%
  \ifthenelse{\equal{\@journal}{Games}}{2073-4336}{}%
  \ifthenelse{\equal{\@journal}{Genes}}{2073-4425}{}%
  \ifthenelse{\equal{\@journal}{Information}}{2078-2489}{}%
  \ifthenelse{\equal{\@journal}{Ijerph}}{1660-4601}{}%
  \ifthenelse{\equal{\@journal}{Ijms}}{1422-0067}{}%
  \ifthenelse{\equal{\@journal}{Jlpea}}{2079-9268}{}%
  \ifthenelse{\equal{\@journal}{MarineDrugs}}{1660-3397}{}%
  \ifthenelse{\equal{\@journal}{Materials}}{1996-1944}{}%
  \ifthenelse{\equal{\@journal}{Micromachines}}{2072-666X}{}%
  \ifthenelse{\equal{\@journal}{Molbank}}{1422-8599}{}%
  \ifthenelse{\equal{\@journal}{Molecules}}{1420-3049}{}%
  \ifthenelse{\equal{\@journal}{Nutrients}}{2072-6643}{}%
  \ifthenelse{\equal{\@journal}{Pharmaceuticals}}{1424-8247}{}%
  \ifthenelse{\equal{\@journal}{Pharmaceutics}}{1999-4923}{}%
  \ifthenelse{\equal{\@journal}{Polymers}}{2073-4360}{}%
  \ifthenelse{\equal{\@journal}{Religions}}{2077-1444}{}%
  \ifthenelse{\equal{\@journal}{RemoteSensing}}{2072-4292}{}%
  \ifthenelse{\equal{\@journal}{Sensors}}{1424-8220}{}%
  \ifthenelse{\equal{\@journal}{Sustainability}}{2071-1050}{}%
  \ifthenelse{\equal{\@journal}{Symmetry}}{2073-8994}{}%
  \ifthenelse{\equal{\@journal}{Toxins}}{2072-6651}{}%
  \ifthenelse{\equal{\@journal}{Viruses}}{1999-4915}{}%
  \ifthenelse{\equal{\@journal}{Water}}{2073-4441}{}%
  }
% % %www.mdpi.com/journal/\logo
\newcommand{\web}{%
  \ifthenelse{\equal{\@journal}{AdministrativeSciences}}{www.mdpi.com/journal/admsci}{}%
  \ifthenelse{\equal{\@journal}{Algorithms}}{www.mdpi.com/journal/algorithms}{}%
  \ifthenelse{\equal{\@journal}{Atmosphere}}{www.mdpi.com/journal/atmosphere}{}%
  \ifthenelse{\equal{\@journal}{Cancers}}{www.mdpi.com/journal/cancers}{}%
  \ifthenelse{\equal{\@journal}{Challenges}}{www.mdpi.com/journal/challenges}{}%
  \ifthenelse{\equal{\@journal}{Computers}}{www.mdpi.com/journal/computers}{}%
  \ifthenelse{\equal{\@journal}{Diagnostics}}{www.mdpi.com/journal/diagnostics}{}%
  \ifthenelse{\equal{\@journal}{Diversity}}{www.mdpi.com/journal/diversity}{}%
  \ifthenelse{\equal{\@journal}{Energies}}{www.mdpi.com/journal/energies}{}%
  \ifthenelse{\equal{\@journal}{Entropy}}{www.mdpi.com/journal/entropy}{}%
  \ifthenelse{\equal{\@journal}{Forests}}{www.mdpi.com/journal/forests}{}%
  \ifthenelse{\equal{\@journal}{FutureInternet}}{www.mdpi.com/journal/futureinternet}{}%
  \ifthenelse{\equal{\@journal}{Games}}{www.mdpi.com/journal/games}{}%
  \ifthenelse{\equal{\@journal}{Genes}}{www.mdpi.com/journal/genes}{}%
  \ifthenelse{\equal{\@journal}{Information}}{www.mdpi.com/journal/information}{}%
  \ifthenelse{\equal{\@journal}{Ijerph}}{www.mdpi.com/journal/ijerph}{}%
  \ifthenelse{\equal{\@journal}{Ijms}}{www.mdpi.com/journal/ijms}{}%
  \ifthenelse{\equal{\@journal}{Jlpea}}{www.mdpi.com/journal/jlpea}{}%
  \ifthenelse{\equal{\@journal}{MarineDrugs}}{www.mdpi.com/journal/marinedrugs}{}%
  \ifthenelse{\equal{\@journal}{Materials}}{www.mdpi.com/journal/materials}{}%
  \ifthenelse{\equal{\@journal}{Micromachines}}{www.mdpi.com/journal/micromachines}{}%
  \ifthenelse{\equal{\@journal}{Molbank}}{www.mdpi.com/journal/molbank}{}%
  \ifthenelse{\equal{\@journal}{Molecules}}{www.mdpi.com/journal/molecules}{}%
  \ifthenelse{\equal{\@journal}{Nutrients}}{www.mdpi.com/journal/nutrients}{}%
  \ifthenelse{\equal{\@journal}{Pharmaceuticals}}{www.mdpi.com/journal/pharmaceuticals}{}%
  \ifthenelse{\equal{\@journal}{Pharmaceutics}}{www.mdpi.com/journal/pharmaceutics}{}%
  \ifthenelse{\equal{\@journal}{Polymers}}{www.mdpi.com/journal/polymers}{}%
  \ifthenelse{\equal{\@journal}{Religions}}{www.mdpi.com/journal/religions}{}%
  \ifthenelse{\equal{\@journal}{RemoteSensing}}{www.mdpi.com/journal/remotesensing}{}%
  \ifthenelse{\equal{\@journal}{Sensors}}{www.mdpi.com/journal/sensors}{}%
  \ifthenelse{\equal{\@journal}{Sustainability}}{www.mdpi.com/journal/sustainability}{}%
  \ifthenelse{\equal{\@journal}{Symmetry}}{www.mdpi.com/journal/symmetry}{}%
  \ifthenelse{\equal{\@journal}{Toxins}}{www.mdpi.com/journal/toxins}{}%
  \ifthenelse{\equal{\@journal}{Viruses}}{www.mdpi.com/journal/viruses}{}%
  \ifthenelse{\equal{\@journal}{Water}}{www.mdpi.com/journal/water}{}%
}
%%%Define the lenght of the header on the front page
\newcommand{\logolength}{%
  \ifthenelse{\equal{\@journal}{AdministrativeSciences}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Algorithms}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Atmosphere}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Cancers}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Challenges}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Computers}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Diagnostics}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Diversity}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Energies}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Entropy}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Forests}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{FutureInternet}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Games}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Genes}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Information}}{\begin{minipage}[t][21mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Ijerph}}{\begin{minipage}[t][25mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Ijms}}{\begin{minipage}[t][21mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Jlpea}}{\begin{minipage}[t][21mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{MarineDrugs}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Materials}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Micromachines}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Molecules}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Molbank}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Nutrients}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Pharmaceuticals}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Pharmaceutics}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Polymers}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Religions}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{RemoteSensing}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Sensors}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Sustainability}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Symmetry}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Toxins}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Viruses}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
  \ifthenelse{\equal{\@journal}{Water}}{\begin{minipage}[t][17mm][t]{0.97\textwidth}}{}%
}
%
%
%Define the header for the front page
\newcommand{\header}{%
\thispagestyle{empty}
\logolength
\begin{flushleft}
  \vspace{-36pt}
  \hspace{-5mm}{\ifthenelse{\equal{\@status}{submit}}{\normalfont Submitted to {\em \journalname}. Pages \thepage \ - \color{black}{\pageref{LastPage}}.
}{\normalfont {\em \journalname}\ {\bfseries {\@pubyear}}, {\em \@pubvolume},
\thepage -\@lastpage; doi:\@doinum}}
\end{flushleft}
%
% Creates command that produces the logo on the first page of an accepted paper.
% For drafts the logo is made in white so that it cannot be seen,
% but leave the formatting unchanged
\vspace{-28pt}
\begin{flushright}
      {\fontsize{9}{9}\selectfont\bfseries{\colorbox{cyan}{\textcolor{white}{OPEN ACCESS}}}}\\[8pt]
\ifthenelse{\equal{\@status}{submit}}{
        {{\logosize}{\logosize}\selectfont\bfseries\itshape \textcolor{white}{\logo}}\\
        \normalfont\bfseries \textcolor{white}{ISSN \issn }\\[-2pt]
        \normalfont \textcolor{white}{www.mdpi.com/journal/\logo}
        }{{\logosize\selectfont\bfseries\itshape \textcolor{logocolor}{\logo}}\\
        \normalfont\bfseries ISSN \issn\\[-2pt]
        \normalfont \web
        }
\end{flushright}
 \end{minipage}}
%
%
%% Maketitle
\renewcommand{\@maketitle}{%
\header
%Title,author, address, abstract, keywords
\begin{flushleft}%
 {%
 \fontsize{12}{12}\selectfont
  \vspace{18pt}
 \textit{\@arttype\\}
 \vspace{6pt}
 }%
  {%3
   \fontsize{18}{18}\selectfont
   \bfseries{\@Title}
   \par
  }%3
  \vspace{12pt}%
  {%4
   \bfseries{\@Author}
   \par
  }%4
  \vspace{12pt}%
{%6
\ifthenelse{\equal{\firstargument}{1}}{}{\hangafter=1\hangindent=0.7em\noindent}%
\@address
\par
}%6
\vspace{12pt}%
{%8
\ifthenelse{\equal{\@authornum}{authors}}{\hangafter=1\hangindent=0.7em\noindent
$^{\star}$ \ifthenelse{\equal{\corresfirstargument}{2}}{Authors}{Author} to whom correspondence should be addressed; \@corres
\par
\vspace{12pt}}{}
}%8
\em
{%9
\ifthenelse{\equal{\@status}{submit}}{Version {\@ \today} submitted to {\em
\@journal}. Typeset by \LaTeX{} using class file mdpi.cls}{\@history}\par
}%9
\end{flushleft}
}%End of maketitle
%%Command for typsetting abstract, kewyords and MSC, PACS
\newcommand{\abstractkeywords}{%
\vspace{12pt}%
\hrule
%
% Specifies how the abstract should be typeset
{%10
\begingroup
 \vspace{6mm}
\leftskip1cm
\rightskip\leftskip
\parindent0em
  \noindent\textbf{Abstract:\space}\@abstract
  \vspace{3mm}
  \par
  \endgroup
}%10
% % Specifies how the keywords should be typeset
% and inserts a horisontal line aterwards
% If PACS and/or MSC are specified they are inserted
% after keywords, but before the line
%
{%11
\begingroup
 \vspace{3pt}
\leftskip1cm
\rightskip\leftskip
\parindent0em
     \textbf{Keywords:\space}\@keyword
 \vspace{12pt}
 \par
 \endgroup
}%11
%
%Specifies how PACS should be typeset
%
{%12
\ifthenelse{\equal{\@pacs}{\@empty}}{}{%
\begingroup
\leftskip1cm
\rightskip\leftskip
\parindent0em
  \noindent\textbf{Classification: PACS\space} \@pacs
  \vspace{12pt}
  \par
  \endgroup}%
}%12
%
% Specifies how PACS should be typeset
%
{%13
\ifthenelse{\equal{\@msc}{\@empty}}{}{%
\begingroup
\leftskip1cm
\rightskip\leftskip
\parindent0em
\noindent\textbf{Classification: MSC\space} \@msc
	\vspace{12pt}
   \par
  \endgroup}%%
}%13
\vspace{6pt}
\hrule
\vspace{27pt}
}%End of Abstract and Keywords
%
%
\AtBeginDocument{%
  \maketitle
  \thispagestyle{empty}
  \let\maketitle\relax
  \ifthenelse{\equal{\@status}{submit}}{\linenumbers}{}
  \abstractkeywords
}
% % % Specifies how section headings should be typeset
\renewcommand{\thesection}{\arabic{section}}
\renewcommand{\section}{ %
  \@startsection{section}{1}{0pt}{12pt}{12pt}{\normalsize\bfseries}}

\renewcommand{\thesubsection}{\arabic{section}.\arabic{subsection}}
\renewcommand{\subsection}{ %
  \@startsection{subsection}{2}{0pt}{12pt}{12pt}{\normalsize\itshape}}

\renewcommand{\thesubsubsection}{}
\renewcommand{\subsubsection}{ %
  \@startsection{subsubsection}{3}{-12pt}{12pt}{12pt}{\normalsize}}

%%%% dot after number of section, but not in ref
\let\@@seccntformat\@seccntformat
\renewcommand*{\@seccntformat}[1]{%
  \expandafter\ifx\csname @seccntformat@#1\endcsname\relax
    \expandafter\@@seccntformat
  \else
    \expandafter
      \csname @seccntformat@#1\expandafter\endcsname
  \fi
    {#1}%
}

\newcommand*{\@seccntformat@subsubsection}[1]{%
  \csname the#1\endcsname\quad
}


\renewcommand*{\@@seccntformat}[1]{%
  \csname the#1\endcsname.\quad
}

%\renewcommand*{\p@section}{\arabic{section}.}
%\renewcommand*{\p@subsection}{\p@section\thesection~}
%\renewcommand*{\p@subsubsection}{\p@subsection\thesubsection~}

%Defines the appendix
\renewcommand{\appendix}{\par
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \setcounter{subsubsection}{0}%
  \gdef\thesection{\@Alph\c@section}%
  \gdef\thesubsection{\@Alph\c@section.\@arabic\c@subsection}%
  \gdef\theHsection{\@Alph\c@section.}%
  \gdef\theHsubsection{\@Alph\c@section.\@arabic\c@subsection}%
  \csname appendixmore\endcsname
}


% Specifies how captions should be typeset
\renewcommand{\@makecaption}[2]{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{\normalsize{\bfseries #1}. #2}%
  \ifdim \wd\@tempboxa >\hsize
    \normalsize{\bfseries #1}. #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}


\setlength{\oddsidemargin}{17.5mm}
\setlength{\voffset}{-1in}
\setlength{\hoffset}{-1in}
\setlength{\topmargin}{17.5mm}
\setlength{\headsep}{0mm}
\setlength{\textheight}{250mm}
\setlength{\textwidth}{175mm}
\setlength{\textfloatsep}{12pt}
\setlength{\abovedisplayskip}{0pt}
\setlength{\belowdisplayskip}{0pt}
\setlength{\parindent}{5mm}
\setlength{\abovecaptionskip}{12pt}
\setlength{\belowcaptionskip}{12pt}
\setlength{\itemsep}{0pt}
\setlength{\partopsep}{0pt}
\setlength{\parsep}{0pt}


% Create headings for pages except the first
\newcommand{\heading}[1]{\gdef\@heading{#1}}
\setlength{\headheight}{8.5pt}
\setlength{\headsep}{25.5pt}

\def\ps@myheadings{%
    \def\@oddfoot{}
    \let\@evenfoot\@oddfoot
    \def\@evenhead{%
     {\ifthenelse{\equal{\@status}{submit}}{Version {\@ \today} submitted to {\em \journalname}%
        \hfil {\bfseries \thepage \ of \pageref{LastPage}}%
}{{\em \journalname} {\bfseries \@pubyear}, {\em \@pubvolume}%
        \hfill {\bfseries \thepage}}%
        \vspace{50pt}}%
       }
    \let\@oddhead\@evenhead
    \let\@mkboth\@gobbletwo
    \let\sectionmark\@gobble
    \let\subsectionmark\@gobble
}

\pagestyle{myheadings}

% Bibliography
\renewcommand\@biblabel[1]{#1. }

\def\thebibliography#1{\vspace{12pt}
\noindent{{\bfseries References}}\vspace{12pt}%
\list{{\arabic{enumi}}}{\def\makelabel##1{\hss{##1}}\topsep=0\p@\parsep=0\p@
\partopsep=0\p@\itemsep=0\p@
\labelsep=1ex
\itemindent=0\p@
\settowidth\labelwidth{\footnotesize[#1]}%
\leftmargin\labelwidth
\advance\leftmargin\labelsep
%\advance\leftmargin -\itemindent
\usecounter{enumi}}
%\def\newblock{\ }
%\sloppy\clubpenalty4000\widowpenalty4000
%\sfcode`\.=1000\relax
}
\let\endthebibliography=\endlist



% Creates command that produces the copyright information.
\newcommand{\cright}{\vspace{12pt}\noindent \copyright
\ifthenelse{\equal{\@status}{submit}}{\ {\@ \today} by the
\@authornum; submitted to {\em \journalname} for possible open
access publication under the terms and conditions of the
Creative Commons Attribution license
http://creativecommons.org/licenses/by/3.0/. }{\ {\@pubyear} by
the \@authornum; licensee MDPI, Basel, Switzerland. This article is an
open access article \linebreak distributed under the terms and conditions
of the Creative Commons Attribution license \linebreak
(http://creativecommons.org/licenses/by/3.0/.)}}

\endinput
