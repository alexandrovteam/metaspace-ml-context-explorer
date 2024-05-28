helm repo add  emblshiny https://sourcecode.embl.de/api/v4/projects/280/packages/helm/stable

helm repo update

helm upgrade --install --namespace metaspace-context-explorer-shinyapp metaspace-context-explorer emblshiny/shiny -f ./values.yaml
