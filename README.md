# Projeto MC859 - Análise de redes de similaridade genética entre doenças: identificação de padrões e comunidades a partir do DisGeNET

## Introdução
  A partir de grandes bases de dados biomédicas é possível modelar relações complexas entre doenças e seus determinantes moleculares como grafos e descobrir padrões não triviais. Neste projeto, utilizaremos o conjunto de dados disponível no repositório [dhimmel/disgenet](https://github.com/dhimmel/disgenet) do GitHub, que contém associações gene-doença processadas a partir do [DisGeNET](https://disgenet.com/).
  
## Objetivos
  O objetivo é identificar clusters de doenças que compartilham mecanismos genéticos, quantificar similaridades entre doenças e produzir visualizações (algo similar a um heatmap nos nós) que evidenciem quais doenças estão mais intimamente relacionadas. Mais especificamente, os objetivos seriam:

- Construir um grafo robusto a partir dos dados processados do repositório (vértices = doenças; arestas = similaridade genética entre doenças).

- Detectar comunidades/_clusters_ de doenças usando métodos de detecção de comunidades.

- Encontrar que tipos de doenças fazem pontes entre cluster.

- Extrair características topológicas dos nós (grau, centralidade, coeficiente de clusterização, betweenness, assortatividade).

- Gerar uma visualização interativa da rede e um heatmap nos nós (cores indicando grau de similaridade/ intensidade da relação).

- Aplicação de _link prediction_ (apagar algumas arestas e tentar adivinhar quais seriam essas arestas apagadas),

- Detectar a presença de _triadic closures_ (a tendência de dois vértices que têm um vizinho em comum se conectarem também entre si) e tentar interpretar esses resultados.

## Dataset

  O projeto utilizará os arquivos disponibilizados no repositório [dhimmel/disgenet](https://github.com/dhimmel/disgenet), que contém associações gene–doença extraídas e processadas a partir do DisGeNET. Esses arquivos já estão organizados em uma tabela, contendo identificadores padronizados para doenças e genes, além de escores de evidência. A escolha por esse repositório se deve ao fato de ele fornecer dados prontos para análise, dispensando a necessidade de implementar consultas diretas à API oficial.
